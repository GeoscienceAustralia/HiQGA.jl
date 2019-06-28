module MCMC_Driver
using TransD_GP, Distributed, SharedArrays, DistributedArrays,
     PyPlot, LinearAlgebra

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
    abs(Temp-1.0) < 1e-12 && TransD_GP.write_history(isample, opt, m, current_misfit[1], stat, wp)
end

function filldarray(a::AbstractArray)
    r = Array{Future, 1}(undef, nworkers())
    for (ipid, pid) in enumerate(workers())
        r[ipid] = @spawnat pid a
    end
    DArray(r)
end

function close_history(wp::DArray)
    for (idx, pid) in enumerate(procs(wp))
        @spawnat pid TransD_GP.close_history(wp[idx])
    end
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

    current_misfit           = DArray(current_misfit_)
    m, opt, stat, opt_EM, d, 
    current_misfit, wp       = map(x -> DArray(x), (m_, opt_, stat_, opt_EM_, d_in_, current_misfit_, wp_))
end

function main(opt_in::TransD_GP.Options, din::AbstractArray, Tmax::Float64, nsamples::Int, opt_EM_in::EMoptions)
    T = 10.0.^range(0, stop = log10(Tmax), length = nworkers())

    till = TransD_GP.closestmultbelow(nsamples-1, opt_in.save_freq) + 1
    nstore = length(1:opt_in.save_freq:till)
    misfit = zeros(Float64, nstore)
    T0store = zeros(Int, nstore)

    m, opt, stat, opt_EM, d, current_misfit, wp = init_chain_darrays(opt_in, opt_EM_in, din[:])

    @show typeof(m)
    @show typeof(d)
    @show typeof(d[1])

    storecount = 1
    t1 = time()
    t2 = time()
    for isample = 1:nsamples

        for ichain in nworkers():-1:2
            jchain = rand(1:ichain)
            if ichain != jchain
                logalpha = (current_misfit[ichain][1] - current_misfit[jchain][1]) *
                                (1.0/T[ichain] - 1.0/T[jchain])
                if log(rand()) < logalpha
                    T[ichain], T[jchain] = T[jchain], T[ichain]
                end
            end
        end

        @sync for(idx, pid) in enumerate(workers())
            @spawnat pid do_mcmc_step(m[idx], opt[idx], stat[idx],
                                    current_misfit[idx], localpart(d),
                                    T[idx], isample, opt_EM[idx], wp[idx])
        end

        if mod(isample-1, 1000) == 0
            dt = time() - t2 #seconds
            t2 = time()
            @info("*****$dt**sec*****")
        end
        if mod(isample-1, opt_in.save_freq) == 0
            T0idx = argmin(abs.(T.-1.0))
            if time() - t1 >= 2. #seconds
                @info("sample: $isample target worker: $(workers()[T0idx]) misfit $(current_misfit[T0idx]) points $(m[T0idx].n)")
                t1 = time()
            end
            #TransD_GP.write_history(isample, opt[T0idx], m[T0idx], current_misfit[T0idx][1], stat[T0idx], wp)
            misfit[storecount] = current_misfit[T0idx][1]
            T0store[storecount] = T0idx
            storecount = storecount + 1
        end
    end

    TransD_GP.close_history(wp)

    return misfit, T0store
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
