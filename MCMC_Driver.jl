module MCMC_Driver
using TransD_GP, Distributed, SharedArrays, DistributedArrays,
     PyPlot, LinearAlgebra, Formatting

mutable struct EMoptions
    sd      :: Float64
    MLnoise :: Bool
end

function EMoptions(;
            sd = 0.0,
            MLnoise = true)
    @assert sd != 0.0
    EMoptions(sd, MLnoise)
end

struct Tpointer
    fp   :: IOStream
    fstr :: String
end

function get_misfit(m::TransD_GP.Model, d::AbstractArray, opt::TransD_GP.Options, opt_EM::EMoptions)
    chi2by2 = 0.0
    select = .!isnan.(d[:])
    if !opt.debug
        r = m.fstar[select] - d[select]
        if !opt_EM.MLnoise
            chi2by2 = 0.5*norm(r/opt_EM.sd)^2
        else
            N = sum(select)
            chi2by2 = 0.5*N*log(norm(r)^2)
        end
    end
    return chi2by2
end

function mh_step!(m::TransD_GP.Model, d::AbstractArray,
    opt::TransD_GP.Options, stat::TransD_GP.Stats,
    Temp::Float64, movetype::Int, current_misfit::Array{Float64, 1}, opt_EM::EMoptions)

    new_misfit = get_misfit(m, d, opt, opt_EM)
    logalpha = (current_misfit[1] - new_misfit)/Temp
    if log(rand()) < logalpha
        current_misfit[1] = new_misfit
        stat.accepted_moves[movetype] += 1
    else
        TransD_GP.undo_move!(movetype, m, opt)
    end
end

function do_mcmc_step(m::TransD_GP.Model, opt::TransD_GP.Options, stat::TransD_GP.Stats,
    current_misfit::Array{Float64, 1}, d::AbstractArray,
    Temp::Float64, isample::Int, opt_EM::EMoptions, wp::TransD_GP.Writepointers)

    # select move and do it
    movetype, priorviolate = TransD_GP.do_move!(m, opt, stat)

    if !priorviolate
        mh_step!(m, d, opt, stat, Temp, movetype, current_misfit, opt_EM)
    end

    # acceptance stats
    TransD_GP.get_acceptance_stats!(isample, opt, stat)

    # write models
    writemodel = false
    abs(Temp-1.0) < 1e-12 && (writemodel = true)
    TransD_GP.write_history(isample, opt, m, current_misfit[1], stat, wp, writemodel)

    return current_misfit[1]
end

function do_mcmc_step(m::DArray{TransD_GP.Model}, opt::DArray{TransD_GP.Options},
    stat::DArray{TransD_GP.Stats}, current_misfit::DArray{Array{Float64, 1}},
    d::AbstractArray, Temp::Float64, isample::Int, opt_EM::DArray{EMoptions},
    wp::DArray{TransD_GP.Writepointers})

    misfit = do_mcmc_step(localpart(m)[1], localpart(opt)[1], localpart(stat)[1],
                            localpart(current_misfit)[1], localpart(d),
                            Temp, isample, localpart(opt_EM)[1], localpart(wp)[1])

end

function close_history(wp::DArray)
    for (idx, pid) in enumerate(procs(wp))
        @spawnat pid TransD_GP.close_history(wp[idx])
    end
end

function open_temperature_file(opt_in::TransD_GP.Options, T::Array{Float64, 1})
    fdataname = opt_in.costs_filename[9:end-4]
    fp_temps  = open(fdataname*"_temps.txt", opt_in.history_mode)
    fmt = "{:d} "
    for i = 1:length(T)-1
        fmt = fmt*"{:f} "
    end
    fmt = fmt*"{:f}"
    tpointer = Tpointer(fp_temps, fmt)
end

function write_temperatures(iter::Int, T::Array{Float64, 1}, tpointer::Tpointer, opt_in::TransD_GP.Options)
    if (mod(iter-1, opt_in.save_freq) == 0 || iter == 1)
        printfmtln(tpointer.fp, tpointer.fstr, iter, T...)
        flush(tpointer.fp)
    end
end

function close_temperature_file(fp::IOStream)
    close(fp)
end

function init_chain_darrays(opt_in::TransD_GP.Options, opt_EM_in::EMoptions, d_in::AbstractArray)
    m_, opt_, stat_, opt_EM_, d_in_, current_misfit_, wp_  = map(x -> Array{Future, 1}(undef, nworkers()), 1:7)

    costs_filename = opt_in.costs_filename[1:end-4]
    fstar_filename = opt_in.fstar_filename[1:end-4]
    x_ftrain_filename = opt_in.x_ftrain_filename[1:end-4]

    for(idx, pid) in enumerate(workers())
        m_[idx]              = @spawnat pid [TransD_GP.init(opt_in)]

        opt_in.costs_filename    = costs_filename*"_$idx.bin"
        opt_in.fstar_filename    = fstar_filename*"_$idx.bin"
        opt_in.x_ftrain_filename = x_ftrain_filename*"_$idx.bin"

        opt_[idx]            = @spawnat pid [opt_in]
        stat_[idx]           = @spawnat pid [TransD_GP.Stats()]
        opt_EM_[idx]         = @spawnat pid [opt_EM_in]
        d_in_[idx]           = @spawnat pid d_in
        current_misfit_[idx] = @spawnat pid [[get_misfit(fetch(m_[idx])[1],
                                              localpart(fetch(d_in_[idx])),
                                              fetch(opt_[idx])[1],
                                              fetch(opt_EM_[idx])[1])]]
        wp_[idx]             = @spawnat pid [TransD_GP.open_history(opt_in)]

    end

    m, opt, stat, opt_EM, d,
    current_misfit, wp       = map(x -> DArray(x), (m_, opt_, stat_, opt_EM_, d_in_, current_misfit_, wp_))
end

function main(opt_in::TransD_GP.Options, din::AbstractArray, Tmax::Float64, nsamples::Int, opt_EM_in::EMoptions)
    T = 10.0.^range(0, stop = log10(Tmax), length = nworkers())

    misfit = zeros(Float64, length(T))

    fp_temps = open_temperature_file(opt_in::TransD_GP.Options, T)
    m, opt, stat, opt_EM, d, current_misfit, wp = init_chain_darrays(opt_in, opt_EM_in, din[:])

    t2 = time()
    for isample = 1:nsamples

        for ichain in nworkers():-1:2
            jchain = rand(1:ichain)
            if ichain != jchain
                logalpha = (misfit[ichain] - misfit[jchain]) *
                                (1.0/T[ichain] - 1.0/T[jchain])
                if log(rand()) < logalpha
                    T[ichain], T[jchain] = T[jchain], T[ichain]
                end
            end
        end

        @sync for(idx, pid) in enumerate(workers())
            @async misfit[idx] = remotecall_fetch(do_mcmc_step, pid, m, opt, stat,
                                    current_misfit, d,
                                    T[idx], isample, opt_EM, wp)
        end

        write_temperatures(isample, T, fp_temps, opt_in)

        if mod(isample-1, 1000) == 0
            dt = time() - t2 #seconds
            t2 = time()
            @info("*****$dt**sec*****")
        end

    end

    close_history(wp)
    close_temperature_file(fp_temps.fp)

    nothing
end

function nicenup(g::PyPlot.Figure;fsize=14)
    for ax in gcf()[:axes]
        ax[:tick_params]("both",labelsize=fsize)
        ax[:xaxis][:label][:set_fontsize](fsize)
        ax[:yaxis][:label][:set_fontsize](fsize)
        ax[:title][:set_fontsize](fsize)
        if typeof(ax[:get_legend_handles_labels]()[1]) != Array{Any,1}
            ax[:legend](loc="best", fontsize=fsize)
        end
    end
    g[:tight_layout]()
end

end
