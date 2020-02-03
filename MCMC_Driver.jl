module MCMC_Driver
using TransD_GP, Distributed, DistributedArrays,
     PyPlot, LinearAlgebra, Formatting

struct Tpointer
    fp   :: IOStream
    fstr :: String
end

function get_misfit(m::TransD_GP.Model, d::AbstractArray, opt::TransD_GP.Options)
    chi2by2 = 0.0
    if !opt.debug
        select = .!isnan.(d[:])
        r = m.fstar[select] - d[select]
        N = sum(select)
        chi2by2 = 0.5*N*log(norm(r)^2)
    end
    return chi2by2
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
    Temp::Float64, movetype::Int, current_misfit::Array{Float64, 1})

    if opt.quasimultid
        new_misfit = get_misfit(m, d, opt, movetype)
    else
        new_misfit = get_misfit(m, d, opt)
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
    Temp::Float64, isample::Int, wp::TransD_GP.Writepointers)

    # select move and do it
    movetype, priorviolate = TransD_GP.do_move!(m, opt, stat)

    if !priorviolate
        mh_step!(m, d, opt, stat, Temp, movetype, current_misfit)
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
    d::AbstractArray, Temp::Float64, isample::Int,
    wp::DArray{TransD_GP.Writepointers})

    misfit = do_mcmc_step(localpart(m)[1], localpart(opt)[1], localpart(stat)[1],
                            localpart(current_misfit)[1], localpart(d),
                            Temp, isample, localpart(wp)[1])

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

function init_chain_darrays(opt_in::TransD_GP.Options, d_in::AbstractArray, chains::Array{Chain, 1})
    m_, opt_, stat_, d_in_, current_misfit_, wp_  = map(x -> Array{Future, 1}(undef, length(chains)), 1:6)

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
        d_in_[idx]           = @spawnat chain.pid d_in
        current_misfit_[idx] = @spawnat chain.pid [[ get_misfit(fetch(m_[idx])[1],
                                               localpart(fetch(d_in_[idx])),
                                               fetch(opt_[idx])[1]) ]]

        wp_[idx]             = @spawnat chain.pid [TransD_GP.open_history(opt_in)]

    end

    m, opt, stat, d,
    current_misfit, wp       = map(x -> DArray(x), (m_, opt_, stat_, d_in_, current_misfit_, wp_))
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

function main(opt_in::TransD_GP.Options, din::AbstractArray;
              nsamples     = 4001,
              nchains      = 1,
              nchainsatone = 1,
              Tmax         = 2.5)


    chains = Chain(nchains, Tmax=Tmax, nchainsatone=nchainsatone)
    m, opt, stat, d, current_misfit, wp = init_chain_darrays(opt_in, din[:], chains)

    t2 = time()
    for isample = 1:nsamples

        swap_temps(chains)

        @sync for(ichain, chain) in enumerate(chains)
            @async chain.misfit = remotecall_fetch(do_mcmc_step, chain.pid, m, opt, stat,
                                                             current_misfit, d,
                                                             chain.T, isample, wp)
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

function getchi2forall(opt_in::TransD_GP.Options;
                        nchains          = 1,
                        figsize          = (17,8),
                        fsize            = 14,
                      )
    if nchains == 1 # then actually find out how many chains there are saved
        nchains = length(filter( x -> occursin(r"misfits.*bin", x), readdir(pwd()) )) # my terrible regex
    end
    # now look at any chain to get how many iterations
    costs_filename = "misfits_"*opt_in.fdataname
    opt_in.costs_filename    = costs_filename*"_1.bin"
    iters          = TransD_GP.history(opt_in, stat=:iter)
    niters         = length(iters)
    # then create arrays of unsorted by temperature T, k, and chi2
    Tacrosschains  = zeros(Float64, niters, nchains)
    kacrosschains  = zeros(Int, niters, nchains)
    X2by2inchains  = zeros(Float64, niters, nchains)
    # get the values into the arrays
    for ichain in 1:nchains
        opt_in.costs_filename = costs_filename*"_$ichain.bin"
        Tacrosschains[:,ichain] = TransD_GP.history(opt_in, stat=:T)
        kacrosschains[:,ichain] = TransD_GP.history(opt_in, stat=:nodes)
        X2by2inchains[:,ichain] = TransD_GP.history(opt_in, stat=:U)
    end

    f, ax = plt.subplots(3,2, sharex=true, figsize=figsize)
    ax[1].plot(iters, kacrosschains)
    ax[1].set_title("unsorted by temperature")
    ax[1].grid()
    ax[1].set_ylabel("# nodes")
    ax[2].plot(iters, X2by2inchains)
    ax[2].grid()
    ax[2].set_ylabel("-Log L")
    ax[3].grid()
    ax[3].plot(iters, Tacrosschains)
    ax[3].set_ylabel("Temperature")
    ax[3].set_xlabel("iterations")

    for jstep = 1:niters
        sortidx = sortperm(vec(Tacrosschains[jstep,:]))
        X2by2inchains[jstep,:] = X2by2inchains[jstep,sortidx]
        kacrosschains[jstep,:] = kacrosschains[jstep,sortidx]
        Tacrosschains[jstep,:] = Tacrosschains[jstep,sortidx]
    end

    nchainsatone = sum(Tacrosschains[1,:] .== 1)
    ax[4].plot(iters, kacrosschains)
    ax[4].set_title("sorted by temperature")
    ax[4].plot(iters, kacrosschains[:,1:nchainsatone], "k")
    ax[4].grid()
    ax[5].plot(iters, X2by2inchains)
    ax[5].plot(iters, X2by2inchains[:,1:nchainsatone], "k")
    ax[5].grid()
    ax[6].plot(iters, Tacrosschains)
    ax[6].plot(iters, Tacrosschains[:,1:nchainsatone], "k")
    ax[6].grid()
    ax[6].set_xlabel("iterations")

    nicenup(f, fsize=fsize)

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
