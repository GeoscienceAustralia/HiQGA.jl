module MCMC_Driver
using TransD_GP, Distributed, DistributedArrays,
     PyPlot, LinearAlgebra, Formatting, UseGA_AEM

mutable struct Sounding
    data :: Array{Float64} 
    x    :: Array{Float64,1}
end

mutable struct EMoptions
    sd        :: Float64
    MLnoise   :: Bool
    soundings :: Array{Sounding, 1}
    ncellsz   :: Int
    dz        :: Float64
end

function EMoptions(;
            sd         = 0.0,
            MLnoise    = true,
            soundings  = nothing,
            ncellsz    = 100,
            dz         = 2.0
                  )
    @assert sd != 0.0
    @assert soundings != nothing
    @assert ncellsz > 10
    @assert dz > 0.5
    EMoptions(sd, MLnoise, soundings, ncellsz, dz)
end

struct Tpointer
    fp   :: IOStream
    fstr :: String
end

function get_misfit(m::TransD_GP.Model, r::AbstractArray, opt::TransD_GP.Options,
                    opt_EM::EMoptions, movetype::Int)
    chi2by2 = 0.0
    if !opt.debug
        recompute = falses(length(opt_EM.soundings))
        for (isounding, sounding) in enumerate(opt_EM.soundings)
            # next line for birth, death, property_change
            δ = m.xtrain_focus[1:end-1,:] - sounding.x
            if norm(δ./opt.influenceradius) < 1.0
                recompute[isounding] = true
            end    
            # next line for position_change
            if movetype == 3 && recompute[isounding] == false
                δ =   m.xtrain_old[1:end-1,:] - sounding.x
                if norm(δ./opt.influenceradius) < 1.0
                    recompute[isounding] = true
                end
            end
            # now for the stride bit in xall
            for s in findall(recompute)
                f = forward(m, soundings[s])
                r[s] = soundings[s].data - f # other stuff *needed* here
            end    
        end
        chi2by2 = sum(r)
    end
    return chi2by2
end

function forward(m::TransD_GP.Model, sounding::Sounding)

end    

function get_misfit(m::TransD_GP.Model, d::AbstractArray, opt::TransD_GP.Options,
                    opt_EM::EMoptions)
    chi2by2 = 0.0
end

mutable struct Chain
    pid           :: Int
    npidsperchain :: Int
    T             :: Float64
    misfit        :: Float64
end

function Chain(nchains::Int;
               Tmax          = 2.5,
               nchainsatone  = 1,
              )

    @assert nchains > 1
    @assert Tmax > 1
    @assert mod(nworkers(), nchains) == 0

    npidsperchain = floor(Int, nworkers()/nchains)
    @info "npidsperchain = $npidsperchain"
    T = 10.0.^range(0, stop = log10(Tmax), length = nchains-nchainsatone+1)
    append!(T, ones(nchainsatone-1))
    chains = Array{Chain, 1}(undef, nchains)

    pid_end = 0
    for ichain in 1:nchains
        pid_start      = pid_end + 1
        pid_end        = pid_start + npidsperchain - 1
        pids           = workers()[pid_start:pid_end]

        chains[ichain] = Chain(pids[1], npidsperchain, T[ichain], 0.0)
    end
    chains
end

function mh_step!(m::TransD_GP.Model, d::AbstractArray,
    opt::TransD_GP.Options, stat::TransD_GP.Stats,
    Temp::Float64, movetype::Int, current_misfit::Array{Float64, 1}, opt_EM::EMoptions)

    if opt.quasimultid
        new_misfit = get_misfit(m, d, opt, opt_EM, movetype)
    else
        new_misfit = get_misfit(m, d, opt, opt_EM)
    end
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
    TransD_GP.write_history(isample, opt, m, current_misfit[1], stat, wp, Temp, writemodel)

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
    @sync for (idx, pid) in enumerate(procs(wp))
        @spawnat pid TransD_GP.close_history(wp[idx])
    end
end

function open_temperature_file(opt_in::TransD_GP.Options, nchains::Int)
    fdataname = opt_in.costs_filename[9:end-4]
    fp_temps  = open(fdataname*"_temps.txt", opt_in.history_mode)
    fmt = "{:d} "
    for i = 1:nchains-1
        fmt = fmt*"{:f} "
    end
    fmt = fmt*"{:f}"
    tpointer = Tpointer(fp_temps, fmt)
end

function write_temperatures(iter::Int, chains::Array{Chain, 1}, tpointer::Tpointer, opt_in::TransD_GP.Options)
    if (mod(iter-1, opt_in.save_freq) == 0 || iter == 1)
        printfmtln(tpointer.fp, tpointer.fstr, iter, (getproperty.(chains,:T))...)
        flush(tpointer.fp)
    end
end

function close_temperature_file(fp::IOStream)
    close(fp)
end

function init_chain_darrays(opt_in::TransD_GP.Options, opt_EM_in::EMoptions, d_in::AbstractArray, chains::Array{Chain, 1})
    m_, opt_, stat_, opt_EM_, d_in_, current_misfit_, wp_  = map(x -> Array{Future, 1}(undef, length(chains)), 1:7)

    costs_filename = "misfits_"*opt_in.fdataname
    fstar_filename = "models_"*opt_in.fdataname
    x_ftrain_filename = "points_"*opt_in.fdataname

    @sync for(idx, chain) in enumerate(chains)
        m_[idx]              = @spawnat chain.pid [TransD_GP.init(opt_in)]

        opt_in.costs_filename    = costs_filename*"_$idx.bin"
        opt_in.fstar_filename    = fstar_filename*"_$idx.bin"
        opt_in.x_ftrain_filename = x_ftrain_filename*"_$idx.bin"

        opt_[idx]            = @spawnat chain.pid [opt_in]
        stat_[idx]           = @spawnat chain.pid [TransD_GP.Stats()]
        opt_EM_[idx]         = @spawnat chain.pid [opt_EM_in]
        d_in_[idx]           = @spawnat chain.pid d_in
        current_misfit_[idx] = @spawnat chain.pid [[ get_misfit(fetch(m_[idx])[1],
                                               localpart(fetch(d_in_[idx])),
                                               fetch(opt_[idx])[1],
                                               fetch(opt_EM_[idx])[1]) ]]

        wp_[idx]             = @spawnat chain.pid [TransD_GP.open_history(opt_in)]

    end

    m, opt, stat, opt_EM, d,
    current_misfit, wp       = map(x -> DArray(x), (m_, opt_, stat_, opt_EM_, d_in_, current_misfit_, wp_))
end

function swap_temps(chains::Array{Chain, 1})
    for ichain in length(chains):-1:2
        jchain = rand(1:ichain)
        if ichain != jchain
            logalpha = (chains[ichain].misfit - chains[jchain].misfit) *
                            (1.0/chains[ichain].T - 1.0/chains[jchain].T)
            if log(rand()) < logalpha
                chains[ichain].T, chains[jchain].T = chains[jchain].T, chains[ichain].T
            end
        end
    end

end

function main(opt_in::TransD_GP.Options, din::AbstractArray, opt_EM_in::EMoptions ;
              nsamples     = 4001,
              nchains      = 1,
              nchainsatone = 1,
              Tmax         = 2.5)


    chains = Chain(nchains, Tmax=Tmax, nchainsatone=nchainsatone)
    m, opt, stat, opt_EM, d, current_misfit, wp = init_chain_darrays(opt_in, opt_EM_in, din[:], chains)

    t2 = time()
    for isample = 1:nsamples

        swap_temps(chains)

        @sync for(ichain, chain) in enumerate(chains)
            @async chain.misfit = remotecall_fetch(do_mcmc_step, chain.pid, m, opt, stat,
                                                             current_misfit, d,
                                                             chain.T, isample, opt_EM, wp)
        end

        if mod(isample-1, 1000) == 0
            dt = time() - t2 #seconds
            t2 = time()
            @info("*****$dt**sec*****")
        end

    end

    close_history(wp)
    nothing
end

function nicenup(g::PyPlot.Figure;fsize=14)
    for ax in gcf().axes
        ax.tick_params("both",labelsize=fsize)
        ax.xaxis.label.set_fontsize(fsize)
        ax.yaxis.label.set_fontsize(fsize)
        ax.title.set_fontsize(fsize)
        if typeof(ax.get_legend_handles_labels()[1]) != Array{Any,1}
            ax.legend(loc="best", fontsize=fsize)
        end
    end
    g.tight_layout()
end

end
