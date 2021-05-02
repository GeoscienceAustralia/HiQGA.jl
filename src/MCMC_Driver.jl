using Distributed, DistributedArrays,
     PyPlot, LinearAlgebra, Formatting, Dates

struct Tpointer
    fp   :: IOStream
    fstr :: String
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

function Chain(chainprocs::Array{Int, 1};
               Tmax          = 2.5,
               nchainsatone  = 1)

    @assert Tmax > 1
    nchains = length(chainprocs)
    npidsperchain = floor(Int, length(chainprocs)/nchains)
    @info "npidsperchain = $npidsperchain"
    T = 10.0.^range(0, stop = log10(Tmax), length = nchains-nchainsatone+1)
    append!(T, ones(nchainsatone-1))
    chains = Array{Chain, 1}(undef, nchains)

    pid_end = 0
    for ichain in 1:nchains
        pid_start      = pid_end + 1
        pid_end        = pid_start + npidsperchain - 1
        pids           = chainprocs[pid_start:pid_end]

        chains[ichain] = Chain(pids[1], npidsperchain, T[ichain], 0.0)
    end
    chains
end

function mh_step!(mns::ModelNonstat, m::ModelStat, mn::ModelNuisance,
    F::Operator, optns::OptionsNonstat, stat::Stats,
    Temp::Float64, movetype::Int, current_misfit::Array{Float64, 1})

    if optns.quasimultid
        if opt.updatenuisances
            new_misfit = get_misfit(mns, mn, optns, movetype, F)
        else
            new_misfit = get_misfit(mns, optns, movetype, F)
        end
    else
        if optns.updatenuisances
            new_misfit = get_misfit(mns, mn, optns, F)
        else
            new_misfit = get_misfit(mns, optns, F)
        end
    end
    logalpha = (current_misfit[1] - new_misfit)/Temp
    if log(rand()) < logalpha
        current_misfit[1] = new_misfit
        stat.accepted_moves[movetype] += 1
    else
        undo_move!(movetype, mns, optns, m)
    end
end

function mh_step!(m::ModelStat, mns::ModelNonstat, mn::ModelNuisance,
    F::Operator, opt::OptionsStat, optns::OptionsNonstat,
    stat::Stats, Temp::Float64, movetype::Int, current_misfit::Array{Float64, 1})

    if opt.quasimultid
        if opt.updatenonstat
            if opt.updatenuisances
                new_misfit = get_misfit(mns, mn, opt, movetype, F)
            else
                new_misfit = get_misfit(mns, opt, movetype, F)
            end
        else
            if opt.updatenuisances
                new_misfit = get_misfit(m, mn, opt, movetype, F)
            else
                new_misfit = get_misfit(m, opt, movetype, F)
            end
        end
    else
        if opt.updatenonstat
            if opt.updatenuisances
                new_misfit = get_misfit(mns, mn, opt, F)
            else
                new_misfit = get_misfit(mns, opt, F)
            end
        else
            if opt.updatenuisances
                new_misfit = get_misfit(m, mn, opt, F)
            else
                new_misfit = get_misfit(m, opt, F)
            end
        end
    end
    logalpha = (current_misfit[1] - new_misfit)/Temp
    if log(rand()) < logalpha
        current_misfit[1] = new_misfit
        stat.accepted_moves[movetype] += 1
    else
        undo_move!(movetype, m, opt, mns, optns)
    end
end

#this only gets called for moves on the nuisance chain
function mh_step!(mn::ModelNuisance, m::ModelStat, mns::ModelNonstat, F::Operator,
    optn::OptionsNuisance, statn::Stats,
    Temp::Float64, movetype::Int, current_misfit::Array{Float64,1})
    if optn.updatenonstat
        new_misfit = get_misfit(mns, mn, optn, F)
    else
        new_misfit = get_misfit(m, mn, optn, F)
    end

    logalpha = (current_misfit[1] - new_misfit)/Temp
    if log(rand()) < logalpha
        current_misfit[1] = new_misfit
        statn.accepted_moves[movetype] += 1
    else
        undo_move!(mn)
    end
end

function do_mcmc_step(m::ModelStat, mns::ModelNonstat, mn::ModelNuisance,
    opt::OptionsStat, optns::OptionsNonstat,
    stat::Stats,
    current_misfit::Array{Float64, 1}, F::Operator,
    Temp::Float64, isample::Int, wp::Writepointers)

    # select move and do it
    movetype, priorviolate = do_move!(m, opt, stat, mns, optns)

    if !priorviolate
        mh_step!(m, mns, mn, F, opt, optns, stat, Temp, movetype, current_misfit)
    end

    # acceptance stats
    get_acceptance_stats!(isample, opt, stat)

    # write models
    writemodel = false
    abs(Temp-1.0) < 1e-12 && (writemodel = true)
    write_history(isample, opt, m, current_misfit[1], stat, wp, Temp, writemodel)

    return current_misfit[1]
end

function do_mcmc_step(m::DArray{ModelStat}, mns::DArray{ModelNonstat},
    mn::DArray{ModelNuisance},
    opt::DArray{OptionsStat}, optns::DArray{OptionsNonstat},
    stat::DArray{Stats}, current_misfit::DArray{Array{Float64, 1}},
    F::DArray{x}, Temp::Float64, isample::Int,
    wp::DArray{Writepointers}) where x<:Operator

    misfit = do_mcmc_step(localpart(m)[1], localpart(mns)[1], localpart(mn)[1],
        localpart(opt)[1], localpart(optns)[1],
        localpart(stat)[1], localpart(current_misfit)[1], localpart(F)[1],
                            Temp, isample, localpart(wp)[1])

end

function do_mcmc_step(mns::ModelNonstat, m::ModelStat, mn::ModelNuisance,
    optns::OptionsNonstat, statns::Stats,
    current_misfit::Array{Float64, 1}, F::Operator,
    Temp::Float64, isample::Int, wp::Writepointers)

    # select move and do it
    movetype, priorviolate = do_move!(mns, m, optns, statns)

    if !priorviolate
        mh_step!(mns, m, mn, F, optns, statns, Temp, movetype, current_misfit)
    end

    # acceptance stats
    get_acceptance_stats!(isample, optns, statns)

    # write models
    writemodel = false
    abs(Temp-1.0) < 1e-12 && (writemodel = true)
    write_history(isample, optns, mns, current_misfit[1], statns, wp, Temp, writemodel)

    return current_misfit[1]
end

function do_mcmc_step(mns::DArray{ModelNonstat}, m::DArray{ModelStat},
    mn::DArray{ModelNuisance},
    optns::DArray{OptionsNonstat}, statns::DArray{Stats},
    current_misfit::DArray{Array{Float64, 1}},
    F::DArray{x}, Temp::Float64, isample::Int,
    wpns::DArray{Writepointers}) where x<:Operator

    misfit = do_mcmc_step(localpart(mns)[1], localpart(m)[1], localpart(mn)[1],
                          localpart(optns)[1], localpart(statns)[1],
                          localpart(current_misfit)[1], localpart(F)[1],
                            Temp, isample, localpart(wpns)[1])

end

function do_mcmc_step(mn::ModelNuisance, m::ModelStat, mns::ModelNonstat,
    optn::OptionsNuisance, statn::Stats, current_misfit::Array{Float64,1},
    F::Operator, Temp::Float64, isample::Int, wpn::Writepointers_nuisance)

    movetype = do_move!(mn, optn, statn)

    mh_step!(mn, m, mns, F, optn, statn, Temp, movetype, current_misfit)

    get_acceptance_stats!(isample, optn, statn)

    writemodel = false
    abs(Temp - 1.0) < 1e-12 && (writemodel = true)
    write_history(isample, optn, mn, current_misfit[1], statn, wpn, Temp, writemodel)

    return current_misfit[1]
end

function do_mcmc_step(mn::DArray{ModelNuisance}, m::DArray{ModelStat},
    mns::DArray{ModelNonstat}, optn::DArray{OptionsNuisance}, statn::DArray{Stats},
    current_misfit::DArray{Array{Float64, 1}},
    F::DArray{x}, Temp::Float64, isample::Int,
    wpn::DArray{Writepointers_nuisance}) where x<:Operator

    misfit = do_mcmc_step(localpart(mn)[1], localpart(m)[1], localpart(mns)[1],
            localpart(optn)[1], localpart(statn)[1], localpart(current_misfit)[1],
            localpart(F)[1], localpart(Temp)[1], localpart(isample)[1], localpart(wpn)[1])
end

function close_history(wp::DArray)
    for (idx, pid) in enumerate(procs(wp))
        @sync @spawnat pid close_history(wp[idx])
    end
end

function open_temperature_file(opt_in::Options, nchains::Int)
    fdataname = opt_in.costs_filename[9:end-4]
    fp_temps  = open(fdataname*"_temps.txt", opt_in.history_mode)
    fmt = "{:d} "
    for i = 1:nchains-1
        fmt = fmt*"{:f} "
    end
    fmt = fmt*"{:f}"
    tpointer = Tpointer(fp_temps, fmt)
end

function write_temperatures(iter::Int, chains::Array{Chain, 1}, tpointer::Tpointer, opt_in::Options)
    if (mod(iter-1, opt_in.save_freq) == 0 || iter == 1)
        printfmtln(tpointer.fp, tpointer.fstr, iter, (getproperty.(chains,:T))...)
        flush(tpointer.fp)
    end
end

function close_temperature_file(fp::IOStream)
    close(fp)
end

function init_chain_darrays(opt_in::OptionsStat,
                            optns_in::OptionsNonstat,
                            optn_in::OptionsNuisance,
                            F_in::Operator, chains::Array{Chain, 1})
    m_, mns_, mn_, opt_, optns_, optn_, F_in_, stat_, statns_, statn_, d_in_,
    current_misfit_, wp_, wpns_, wpn_  = map(x -> Array{Future, 1}(undef, length(chains)), 1:15)

    costs_filename = "misfits_"*opt_in.fdataname
    fstar_filename = "models_"*opt_in.fdataname
    x_ftrain_filename = "points_"*opt_in.fdataname
    nu_filename = "values_nuisance_"*opt_in.fdataname

    if opt_in.history_mode == "a"
        setrestartflag.([opt_in, optns_in, optn_in])
    end

    iterlast = 0
    @sync for(idx, chain) in enumerate(chains)

        opt_in.costs_filename      = costs_filename*"s_$idx.bin"
        optns_in.costs_filename    = costs_filename*"ns_$idx.bin"
        opt_in.fstar_filename      = fstar_filename*"s_$idx.bin"
        optns_in.fstar_filename    = fstar_filename*"ns_$idx.bin"
        opt_in.x_ftrain_filename   = x_ftrain_filename*"s_$idx.bin"
        optns_in.x_ftrain_filename = x_ftrain_filename*"ns_$idx.bin"
        optn_in.costs_filename     = costs_filename*"nuisance_$idx.bin"
        optn_in.vals_filename      = nu_filename*"$idx.bin"

        opt_[idx]            = @spawnat chain.pid [opt_in]
        optns_[idx]          = @spawnat chain.pid [optns_in]
        optn_[idx]           = @spawnat chain.pid [optn_in]

        m_[idx]                = @spawnat chain.pid [init(opt_in)]
        mns_[idx]              = @spawnat chain.pid [init(optns_in,
                                                            fetch(m_[idx])[1])]
        mn_[idx]               = @spawnat chain.pid [init(optn_in)]

        @sync wp_[idx]             = @spawnat chain.pid [open_history(opt_in)]
        #@info "sending $(optns_in.costs_filename)"
        @sync wpns_[idx]           = @spawnat chain.pid [open_history(optns_in)]
        @sync wpn_[idx]            = @spawnat chain.pid [open_history(optn_in)]

        stat_[idx]           = @spawnat chain.pid [Stats()]
        statns_[idx]         = @spawnat chain.pid [Stats()]
        statn_[idx]          = @spawnat chain.pid [Stats(nmoves=optn_in.nnu)]

        F_in_[idx]           = @spawnat chain.pid [F_in]
        if opt_in.updatenonstat
            if opt_in.updatenuisances
                current_misfit_[idx] = @spawnat chain.pid [[ get_misfit(fetch(mns_[idx])[1],
                                               fetch(mn_[idx])[1],
                                               fetch(optns_[idx])[1],
                                               fetch(F_in_[idx])[1]) ]]
            else
                current_misfit_[idx] = @spawnat chain.pid [[ get_misfit(fetch(mns_[idx])[1],
                                               fetch(optns_[idx])[1],
                                               fetch(F_in_[idx])[1]) ]]
            end
        else
            if opt_in.updatenuisances
                current_misfit_[idx] = @spawnat chain.pid [[ get_misfit(fetch(m_[idx])[1],
                                               fetch(mn_[idx])[1],
                                               fetch(opt_[idx])[1],
                                               fetch(F_in_[idx])[1]) ]]
            else
                current_misfit_[idx] = @spawnat chain.pid [[ get_misfit(fetch(m_[idx])[1],
                                               fetch(opt_[idx])[1],
                                               fetch(F_in_[idx])[1]) ]]
            end
        end
        if opt_in.history_mode=="a"
            if idx == length(chains)
                iterlast = history(opt_in, stat=:iter)[end]
            end
            chains[idx].T = history(opt_in, stat=:T)[end]
        end
    end
    m, mns, mn, opt, optns, optn, stat, statns, statn, F,
    current_misfit, wp, wpns, wpn = map(x -> DArray(x), (m_, mns_, mn_, opt_, optns_, optn_,
                                    stat_, statns_, statn_, F_in_, current_misfit_,
                                    wp_, wpns_, wpn_))
#return mn as well (at least, possibly return some "nuisance stat" as well)
    @info "initialisation complete"
    return m, mns, mn, opt, optns, optn, stat, statns, statn, F, current_misfit,
            wp, wpns, wpn, iterlast
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

# non pbs: this is for everything nonstationary, stationary and nuisance
# everyday users should not call this but the
# versions underneath instead
function main(opt_in       ::OptionsStat,
              optns_in     ::OptionsNonstat,
              optn_in      ::OptionsNuisance,
              F_in         ::Operator;
              nsamples     = 4001,
              nchains      = 1,
              nchainsatone = 1,
              Tmax         = 2.5)


    chains = Chain(nchains, Tmax=Tmax, nchainsatone=nchainsatone)
    m, mns, mn, opt, optns, optn, stat, statns, statn,
    F, current_misfit, wp, wpns, wpn, iterlast = init_chain_darrays(opt_in,
                                                optns_in, optn_in, F_in, chains)

    domcmciters(iterlast, nsamples, chains, opt_in, mns, m, mn, optns, opt,
                optn, statns, stat, statn, current_misfit, F, wpns, wp, wpn)

    close_history(wp)
    close_history(wpns)
    close_history(wpn)
    d_closeall()
    nothing
end

# non pbs: this is for stationary only
function main(opt_in       ::OptionsStat,
              F_in         ::Operator;
              nsamples     = 4001,
              nchains      = 1,
              nchainsatone = 1,
              Tmax         = 2.5)

    @assert opt_in.needλ²fromlog == false
    @assert opt_in.updatenonstat == false
    updatenuisances = false

    optns_dummy = transD_GP.OptionsNonstat(opt_in,
            nmin = 2,
            nmax = 3,
            fbounds = opt_in.fbounds,
            δ = opt_in.δ,
            demean = opt_in.demean,
            sdev_prop = opt_in.sdev_prop,
            sdev_pos = opt_in.sdev_pos,
            pnorm = opt_in.pnorm,
            K = opt_in.K
            )

    optn_dummy = OptionsNuisance([0.], [0. 0], 0, updatenuisances,
     opt_in.updatenonstat, opt_in.debug, opt_in.stat_window, opt_in.dispstatstoscreen,
     opt_in.report_freq, opt_in.save_freq, opt_in.history_mode, "misfits_nuisance_"*opt_in.fdataname*".bin",
     "values_nuisance_"*opt_in.fdataname*".bin", opt_in.fdataname, [0 0.], [0], [0. 0.], [0.])


    main(opt_in, optns_dummy, optn_dummy, F_in,
            nsamples     = nsamples,
            nchains      = nchains,
            nchainsatone = nchainsatone,
            Tmax         = Tmax)
end

# non pbs: this is for stationary and nuisance only
function main(opt_in       ::OptionsStat,
              optn_in      ::OptionsNuisance,
              F_in         ::Operator;
              nsamples     = 4001,
              nchains      = 1,
              nchainsatone = 1,
              Tmax         = 2.5)

    @assert opt_in.needλ²fromlog == false
    @assert opt_in.updatenonstat == false
    @assert optn_in.updatenuisances == true
    @assert opt_in.updatenuisances == true

    optdummy = transD_GP.OptionsNonstat(opt_in,
            nmin = 2,
            nmax = 3,
            fbounds = opt_in.fbounds,
            δ = opt_in.δ,
            demean = opt_in.demean,
            sdev_prop = opt_in.sdev_prop,
            sdev_pos = opt_in.sdev_pos,
            pnorm = opt_in.pnorm,
            K = opt_in.K
            )

    main(opt_in, optdummy, optn_in, F_in,
            nsamples     = nsamples,
            nchains      = nchains,
            nchainsatone = nchainsatone,
            Tmax         = Tmax)
end

# non pbs: this is for stationary and nonstationary only
function main(opt_in       ::OptionsStat,
              optns_in     ::OptionsNonstat,
              F_in         ::Operator;
              nsamples     = 4001,
              nchains      = 1,
              nchainsatone = 1,
              Tmax         = 2.5)

    @assert opt_in.needλ²fromlog == true
    @assert opt_in.updatenonstat == true
    updatenuisances = false

    optdummy = OptionsNuisance([0.], [0. 0], 0, updatenuisances,
     opt_in.updatenonstat, opt_in.debug, opt_in.stat_window, opt_in.dispstatstoscreen,
     opt_in.report_freq, opt_in.save_freq, opt_in.history_mode, "misfits_nuisance_"*opt_in.fdataname*".bin",
     "values_nuisance_"*opt_in.fdataname*".bin", opt_in.fdataname, [0 0.], [0], [0. 0.], [0.])

    main(opt_in, optns_in, optdummy, F_in,
            nsamples     = nsamples,
            nchains      = nchains,
            nchainsatone = nchainsatone,
            Tmax         = Tmax)
end

# non pbs: If certainly wanting to do stat, nonstat and nuisance all together
function main(doall        ::Bool,
              opt_in       ::OptionsStat,
              optns_in     ::OptionsNonstat,
              optn_in      ::OptionsNuisance,
              F_in         ::Operator;
              nsamples     = 4001,
              nchains      = 1,
              nchainsatone = 1,
              Tmax         = 2.5)
    @assert doall == true
    @assert opt_in.needλ²fromlog == true
    @assert opt_in.updatenonstat == true
    @assert optn_in.updatenuisances == true

    main(opt_in, optns_in, optn_in, F_in,
            nsamples     = nsamples,
            nchains      = nchains,
            nchainsatone = nchainsatone,
            Tmax         = Tmax)
end

# This is for the cluster,
# again ordinary users should not use this but versions underneath
function main(opt_in       ::OptionsStat,
              optns_in     ::OptionsNonstat,
              optn_in      ::OptionsNuisance,
              F_in         ::Operator,
              chainprocs   ::Array{Int, 1};
              nsamples     = 4001,
              nchainsatone = 1,
              Tmax         = 2.5)

    chains = Chain(chainprocs, Tmax=Tmax)
    m, mns, mn, opt, optns, optn, stat, statns, statn,
    F, current_misfit, wp, wpns, wpn, iterlast = init_chain_darrays(opt_in,
                                                optns_in, optn_in, F_in, chains)

    domcmciters(iterlast, nsamples, chains, opt_in, mns, m, mn, optns, opt,
                optn, statns, stat, statn, current_misfit, F, wpns, wp, wpn)

    close_history(wp)
    close_history(wpns)
    close_history(wpn)
    nothing
end

# on cluster: this is for stationary only
function main(opt_in       ::OptionsStat,
              F_in         ::Operator,
              chainprocs   ::Array{Int, 1};
              nsamples     = 4001,
              nchainsatone = 1,
              Tmax         = 2.5)

    @assert opt_in.needλ²fromlog == false
    @assert opt_in.updatenonstat == false
    updatenuisances = false

    optns_dummy = transD_GP.OptionsNonstat(opt_in,
            nmin = 2,
            nmax = 3,
            fbounds = opt_in.fbounds,
            δ = opt_in.δ,
            demean = opt_in.demean,
            sdev_prop = opt_in.sdev_prop,
            sdev_pos = opt_in.sdev_pos,
            pnorm = opt_in.pnorm,
            K = opt_in.K
            )

    optn_dummy = OptionsNuisance([0.], [0. 0], 0, updatenuisances,
     opt_in.updatenonstat, opt_in.debug, opt_in.stat_window, opt_in.dispstatstoscreen,
     opt_in.report_freq, opt_in.save_freq, opt_in.history_mode, "misfits_nuisance_"*opt_in.fdataname*".bin",
     "values_nuisance_"*opt_in.fdataname*".bin", opt_in.fdataname, [0 0.], [0], [0. 0.], [0.])


    main(opt_in, optns_dummy, optn_dummy, F_in, chainprocs,
            nsamples     = nsamples,
            nchainsatone = nchainsatone,
            Tmax         = Tmax)
end

# on cluster: stat + nuisance only
function main(opt_in       ::OptionsStat,
              optn_in      ::OptionsNuisance,
              F_in         ::Operator,
              chainprocs   ::Array{Int, 1};
              nsamples     = 4001,
              nchainsatone = 1,
              Tmax         = 2.5)

    @assert opt_in.needλ²fromlog == false
    @assert opt_in.updatenonstat == false
    @assert optn_in.updatenuisances == true
    @assert opt_in.updatenuisances == true

    optdummy = transD_GP.OptionsNonstat(opt_in,
            nmin = 2,
            nmax = 3,
            fbounds = opt_in.fbounds,
            δ = opt_in.δ,
            demean = opt_in.demean,
            sdev_prop = opt_in.sdev_prop,
            sdev_pos = opt_in.sdev_pos,
            pnorm = opt_in.pnorm,
            K = opt_in.K
            )

    main(opt_in, optdummy, optn_in, F_in, chainprocs,
            nsamples     = nsamples,
            nchainsatone = nchainsatone,
            Tmax         = Tmax)
end

# on cluster: stationary and nonstationary only
function main(opt_in       ::OptionsStat,
              optns_in     ::OptionsNonstat,
              F_in         ::Operator,
              chainprocs   ::Array{Int, 1};
              nsamples     = 4001,
              nchainsatone = 1,
              Tmax         = 2.5)

    @assert opt_in.needλ²fromlog == true
    @assert opt_in.updatenonstat == true
    updatenuisances = false

    optdummy = OptionsNuisance([0.], [0. 0], 0, updatenuisances,
     opt_in.updatenonstat, opt_in.debug, opt_in.stat_window, opt_in.dispstatstoscreen,
     opt_in.report_freq, opt_in.save_freq, opt_in.history_mode, "misfits_nuisance_"*opt_in.fdataname*".bin",
     "values_nuisance_"*opt_in.fdataname*".bin", opt_in.fdataname, [0 0.], [0], [0. 0.], [0.])

    main(opt_in, optns_in, optdummy, F_in, chainprocs,
            nsamples     = nsamples,
            nchainsatone = nchainsatone,
            Tmax         = Tmax)
end

# on cluster: stat, nonstat and nuisance all together
function main(doall        ::Bool,
              opt_in       ::OptionsStat,
              optns_in     ::OptionsNonstat,
              optn_in      ::OptionsNuisance,
              F_in         ::Operator,
              chainprocs   ::Array{Int, 1};
              nsamples     = 4001,
              nchainsatone = 1,
              Tmax         = 2.5)
    @assert doall == true
    @assert opt_in.needλ²fromlog == true
    @assert opt_in.updatenonstat == true
    @assert optn_in.updatenuisances == true

    main(opt_in, optns_in, optn_in, F_in, chainprocs,
            nsamples     = nsamples,
            nchainsatone = nchainsatone,
            Tmax         = Tmax)
end

# all main() call this
function domcmciters(iterlast, nsamples, chains, opt_in, mns, m, mn, optns, opt,
            optn, statns, stat, statn, current_misfit, F, wpns, wp, wpn)
    t2 = time()
    for isample = iterlast+1:iterlast+nsamples

        swap_temps(chains)
        if opt_in.updatenonstat
            @sync for(ichain, chain) in enumerate(chains)
                @async chain.misfit = remotecall_fetch(do_mcmc_step, chain.pid,
                                                mns, m, mn, optns, statns,
                                                current_misfit, F,
                                                chain.T, isample, wpns)
            end
        end
        if opt_in.updatenuisances
            @sync for(ichain, chain) in enumerate(chains)
                @async chain.misfit = remotecall_fetch(do_mcmc_step, chain.pid,
                                                mn, m, mns, optn, statn,
                                                current_misfit, F,
                                                chain.T, isample, wpn)
            end
        end
        @sync for(ichain, chain) in enumerate(chains)
            @async chain.misfit = remotecall_fetch(do_mcmc_step, chain.pid,
                                            m, mns, mn, opt, optns, stat,
                                            current_misfit, F,
                                            chain.T, isample, wp)
        end
        if mod(isample-1, 1000) == 0
            dt = time() - t2 #seconds
            t2 = time()
            @info("*****$dt**sec*****")
        end
    end
end

# maybe split this off into a different include file
# no nuisances e.g., SkyTEM
function loopacrosssoundings(soundings::Array{S, 1}, opt_in::Options;
                            nsequentialiters   =-1,
                            nparallelsoundings =-1,
                            zfixed             = [-1e5],
                            ρfixed             = [1e12],
                            useML              = false,
                            zstart             = 0.0,
                            extendfrac         = 1.06,
                            dz                 = 2.,
                            ρbg                = 10,
                            nlayers            = 50,
                            ntimesperdecade    = 10,
                            nfreqsperdecade    = 5,
                            Tmax               = -1,
                            nsamples           = -1,
                            nchainsatone       = -1,
                            nchainspersounding = -1) where S<:Sounding

    @assert nsequentialiters  != -1
    @assert nparallelsoundings != -1
    @assert nchainspersounding != -1
    @assert nsamples != - 1
    @assert nchainsatone != -1
    @assert Tmax != -1

    nsoundings = length(soundings)
    opt= deepcopy(opt_in)

    for iter = 1:nsequentialiters
        if iter<nsequentialiters
            ss = (iter-1)*nparallelsoundings+1:iter*nparallelsoundings
        else
            ss = (iter-1)*nparallelsoundings+1:nsoundings
        end
        @info "soundings in loop $iter of $nsequentialiters", ss
        r_nothing = Array{Nothing, 1}(undef, length(ss))
        @sync for (i, s) in Iterators.reverse(enumerate(ss))
            pids = (i-1)*nchainspersounding+i:i*(nchainspersounding+1)
            @info "pids in sounding $s:", pids

            aem, = makeoperator(    soundings[s],
                                    zfixed = zfixed,
                                    ρfixed = ρfixed,
                                    zstart = zstart,
                                    extendfrac = extendfrac,
                                    dz = dz,
                                    ρbg = ρbg,
                                    useML = useML,
                                    nlayers = nlayers,
                                    ntimesperdecade = ntimesperdecade,
                                    nfreqsperdecade = nfreqsperdecade)

            opt = deepcopy(opt_in)
            opt.fdataname = soundings[s].sounding_string*"_"

            @async r_nothing[i] = remotecall_fetch(main, pids[1], opt, aem, collect(pids[2:end]),
                                    Tmax         = Tmax,
                                    nsamples     = nsamples,
                                    nchainsatone = nchainsatone)

        end # @sync
        @info "done $iter out of $nsequentialiters at $(Dates.now())"
    end
end

# there are definitely nuisances, e.g. TEMPEST
function loopacrosssoundings(soundings::Array{S, 1};
                                nsequentialiters   = -1,
                                nparallelsoundings = -1,
                                zfixed             = [-1e5],
                                ρfixed             = [1e12],
                                zstart             = 0.0,
                                extendfrac         = 1.06,
                                dz                 = 2.,
                                ρbg                = 10,
                                nlayers            = 50,
                                ntimesperdecade    = 10,
                                nfreqsperdecade    = 5,
                                Tmax               = -1,
                                nsamples           = -1,
                                nchainsatone       = -1,
                                nchainspersounding = -1,
                                znall              = [1],
                                nmin               = 2,
                                nmax               = 40,
                                K                  = GP.Mat32(),
                                demean             = true,
                                sdpos              = 0.05,
                                sdprop             = 0.05,
                                fbounds            = [-0.5 2.5],
                                λ                  = [2],
                                δ                  = 0.1,
                                pnorm              = 2,
                                save_freq          = 50,
                                nuisance_sdev      = [0.],
                                nuisance_bounds    = [0. 0.],
                                updatenuisances    = true,
                                dispstatstoscreen  = false,
                                useML              = false,
                                restart            = false) where S<:Sounding

    @assert nsequentialiters  != -1
    @assert nparallelsoundings != -1
    @assert nchainspersounding != -1
    @assert nsamples != - 1
    @assert nchainsatone != -1
    @assert Tmax != -1

    nsoundings = length(soundings)

    for iter = 1:nsequentialiters
        if iter<nsequentialiters
            ss = (iter-1)*nparallelsoundings+1:iter*nparallelsoundings
        else
            ss = (iter-1)*nparallelsoundings+1:nsoundings
        end
        @info "soundings in loop $iter of $nsequentialiters", ss
        r_nothing = Array{Nothing, 1}(undef, length(ss))
        @sync for (i, s) in Iterators.reverse(enumerate(ss))
            pids = (i-1)*nchainspersounding+i:i*(nchainspersounding+1)
            @info "pids in sounding $s:", pids

            aem, znall = makeoperator(    soundings[s],
                                    zfixed = zfixed,
                                    ρfixed = ρfixed,
                                    zstart = zstart,
                                    extendfrac = extendfrac,
                                    dz = dz,
                                    ρbg = ρbg,
                                    useML = useML,
                                    nlayers = nlayers,
                                    ntimesperdecade = ntimesperdecade,
                                    nfreqsperdecade = nfreqsperdecade)

            opt, optn = transD_GP.make_tdgp_opt(soundings[s],
                                znall = znall,
                                fileprefix = soundings[s].sounding_string,
                                nmin = nmin,
                                nmax = nmax,
                                K = K,
                                demean = demean,
                                sdpos = sdpos,
                                sdprop = sdprop,
                                fbounds = fbounds,
                                save_freq = save_freq,
                                λ = λ,
                                δ = δ,
                                nuisance_bounds = nuisance_bounds,
                                nuisance_sdev = nuisance_sdev,
                                updatenuisances = updatenuisances,
                                restart = restart,
                                dispstatstoscreen = dispstatstoscreen
                                )

            @async r_nothing[i] = remotecall_fetch(main, pids[1], opt, optn, aem, collect(pids[2:end]),
                                    Tmax         = Tmax,
                                    nsamples     = nsamples,
                                    nchainsatone = nchainsatone)

        end # @sync
        @info "done $iter out of $nsequentialiters at $(Dates.now())"
    end
end
