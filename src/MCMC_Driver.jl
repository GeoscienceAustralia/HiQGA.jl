using Distributed, DistributedArrays,
     PyPlot, LinearAlgebra, Formatting, Dates

import .AbstractOperator.Operator
import .AbstractOperator.Sounding
import .AbstractOperator.Operator
import .AbstractOperator.get_misfit
import .AbstractOperator.Sounding
import .AbstractOperator.makeoperator
import .AbstractOperator.make_tdgp_opt
import .AbstractOperator.getresidual

struct Tpointer
    fp   :: IOStream
    fstr :: String
end

mutable struct Chain
    master_pid    :: Int
    pid           :: Int
    npidsperchain :: Int
    T             :: Float64
    misfit        :: Float64
end

function Chain(
               master_pid::Int,
               chainprocs::Array{Int, 1};
               Tmax          = 2.5,
               nchainsatone  = 1)

    @assert Tmax > 1
    nchains = length(chainprocs)
    npidsperchain = floor(Int, length(chainprocs)/nchains)
    T = 10.0.^range(0, stop = log10(Tmax), length = nchains-nchainsatone+1)
    append!(T, ones(nchainsatone-1))
    chains = Array{Chain, 1}(undef, nchains)

    pid_end = 0
    for ichain in 1:nchains
        pid_start      = pid_end + 1
        pid_end        = pid_start + npidsperchain - 1
        pids           = chainprocs[pid_start:pid_end]

        chains[ichain] = Chain(master_pid, pids[1], npidsperchain, T[ichain], 0.0)
    end
    chains
end

function accept(current_misfit, new_misfit, stat, Temp, movetype) 
    logalpha = (current_misfit[1] - new_misfit)/Temp
    accepted = true
    if log(rand()) < logalpha
        current_misfit[1] = new_misfit
        stat.accepted_moves[movetype] += 1
    else
        accepted = false
    end
    accepted
end

function mh_step!(m::ModelStat, mn::ModelNuisance,
    F::Operator, opt::OptionsStat, stat::Stats,
    Temp::Float64, movetype::Int, current_misfit::Array{Float64, 1})
    # for purely stat move and nuisance
    new_misfit = get_misfit(m, mn, opt, F)
    accepted = accept(current_misfit, new_misfit, stat, Temp, movetype)
    !accepted && undo_move!(movetype, m, opt)
end

function mh_step!(mns::ModelNonstat, m::ModelStat, mn::ModelNuisance,
    F::Operator, optns::OptionsNonstat, stat::Stats,
    Temp::Float64, movetype::Int, current_misfit::Array{Float64, 1})
    # for nonstat move and nuisance
    new_misfit = get_misfit(mns, mn, optns, F)
    accepted = accept(current_misfit, new_misfit, stat, Temp, movetype)
    !accepted && undo_move!(movetype, mns, optns, m)
end

function mh_step!(m::ModelStat, mns::ModelNonstat, mn::ModelNuisance,
    F::Operator, opt::OptionsStat, optns::OptionsNonstat,
    stat::Stats, Temp::Float64, movetype::Int, current_misfit::Array{Float64, 1})
    # for stat move updating nonstat and nuisance
    new_misfit = get_misfit(mns, mn, opt, F)
    accepted = accept(current_misfit, new_misfit, stat, Temp, movetype)
    !accepted && undo_move!(movetype, m, opt, mns, optns)
end

function mh_step!(mn::ModelNuisance, m::Model, F::Operator,
    optn::OptionsNuisance, statn::Stats,
    Temp::Float64, movetype::Int, current_misfit::Array{Float64,1})
    # this only gets called for moves on the nuisance chain
    # can be either nonstat or stat
    new_misfit = get_misfit(m, mn, optn, F)
    accepted = accept(current_misfit, new_misfit, statn, Temp, movetype)
    !accepted && undo_move!(mn)
end

function mh_step!(m::ModelStat, F::Operator, opt::OptionsStat, stat::Stats,
    Temp::Float64, movetype::Int, current_misfit::Array{Float64, 1})
    # for purely stat move
    new_misfit = get_misfit(m, opt, F)
    accepted = accept(current_misfit, new_misfit, stat, Temp, movetype)
    !accepted && undo_move!(movetype, m, opt)
end

function mh_step!(mns::ModelNonstat, m::ModelStat, 
    F::Operator, optns::OptionsNonstat, statns::Stats,
    Temp::Float64, movetype::Int, current_misfit::Array{Float64, 1})
    # for purely nonstat move 
    new_misfit = get_misfit(mns, optns, F)
    accepted = accept(current_misfit, new_misfit, statns, Temp, movetype)
    !accepted && undo_move!(movetype, mns, optns, m)
end

function mh_step!(m::ModelStat, mns::ModelNonstat,
    F::Operator, opt::OptionsStat, optns::OptionsNonstat,
    stat::Stats, Temp::Float64, movetype::Int, current_misfit::Array{Float64, 1})
    # for stat move updating nonstat
    new_misfit = get_misfit(mns, optns, F)
    accepted = accept(current_misfit, new_misfit, stat, Temp, movetype)
    !accepted && undo_move!(movetype, m, opt, mns, optns)
end

function do_mcmc_step(m::ModelStat, mns::ModelNonstat, mn::ModelNuisance,
    opt::OptionsStat, optns::OptionsNonstat, stat::Stats,
    current_misfit::Array{Float64, 1}, F::Operator,
    Temp::Float64, isample::Int, wp::Writepointers, chain_idx::Int, master_pid::Int)
    # Stationary GP changes which update nonstationary GP + nuisance
    movetype, priorviolate = do_move!(m, opt, stat, mns, optns)
    if !priorviolate
        mh_step!(m, mns, mn, F, opt, optns, stat, Temp, movetype, current_misfit)
    end
    get_acceptance_stats!(isample, opt, stat)
    writemodel = false
    abs(Temp-1.0) < 1e-12 && (writemodel = true)
    write_history(isample, opt, m, current_misfit[1], stat, wp, Temp, writemodel, chain_idx, master_pid)
    return current_misfit[1]
end

function do_mcmc_step(m::DArray{ModelStat}, mns::DArray{ModelNonstat},
    mn::DArray{ModelNuisance}, opt::DArray{OptionsStat}, 
    optns::DArray{OptionsNonstat}, stat::DArray{Stats}, 
    current_misfit::DArray{Array{Float64, 1}}, F::DArray{x}, Temp::Float64, 
    isample::Int, wp::Writepointers, chain_idx::Int, master_pid::Int) where x<:Operator
    misfit = do_mcmc_step(localpart(m)[1], localpart(mns)[1], localpart(mn)[1],
        localpart(opt)[1], localpart(optns)[1], localpart(stat)[1], 
        localpart(current_misfit)[1], localpart(F)[1],
        Temp, isample, wp, chain_idx, master_pid)
end

function do_mcmc_step(mns::ModelNonstat, m::ModelStat, mn::ModelNuisance,
    optns::OptionsNonstat, statns::Stats, current_misfit::Array{Float64, 1},
    F::Operator, Temp::Float64, isample::Int, wp::Writepointers, chain_idx::Int, master_pid::Int)
    # purely nonstationary GP moves + nuisance
    movetype, priorviolate = do_move!(mns, m, optns, statns)
    if !priorviolate
        mh_step!(mns, m, mn, F, optns, statns, Temp, movetype, current_misfit)
    end
    get_acceptance_stats!(isample, optns, statns)
    writemodel = false
    abs(Temp-1.0) < 1e-12 && (writemodel = true)
    write_history(isample, optns, mns, current_misfit[1], statns, wp, Temp, writemodel, chain_idx, master_pid)
    return current_misfit[1]
end

function do_mcmc_step(mns::DArray{ModelNonstat}, m::DArray{ModelStat},
    mn::DArray{ModelNuisance}, optns::DArray{OptionsNonstat}, statns::DArray{Stats},
    current_misfit::DArray{Array{Float64, 1}}, F::DArray{x}, Temp::Float64, 
    isample::Int, wpns::Writepointers, chain_idx::Int, master_pid::Int) where x<:Operator
    misfit = do_mcmc_step(localpart(mns)[1], localpart(m)[1], localpart(mn)[1],
        localpart(optns)[1], localpart(statns)[1],localpart(current_misfit)[1], 
        localpart(F)[1], Temp, isample, wpns, chain_idx, master_pid)
end

function do_mcmc_step(m::ModelStat, mn::ModelNuisance,
    opt::OptionsStat, stat::Stats, current_misfit::Array{Float64, 1},
    F::Operator, Temp::Float64, isample::Int, wp::Writepointers, chain_idx::Int, master_pid::Int)
    # purely stationary GP moves + nuisance
    movetype, priorviolate = do_move!(m, opt, stat)
    if !priorviolate
        mh_step!(m, mn, F, opt, stat, Temp, movetype, current_misfit)
    end
    get_acceptance_stats!(isample, opt, stat)
    writemodel = false
    abs(Temp-1.0) < 1e-12 && (writemodel = true)
    write_history(isample, opt, m, current_misfit[1], stat, wp, Temp, writemodel, chain_idx, master_pid)
    return current_misfit[1]
end

function do_mcmc_step(m::DArray{ModelStat}, mn::DArray{ModelNuisance}, 
    opt::DArray{OptionsStat}, stat::DArray{Stats},
    current_misfit::DArray{Array{Float64, 1}}, F::DArray{x}, Temp::Float64, 
    isample::Int, wp::Writepointers, chain_idx::Int, master_pid::Int) where x<:Operator
    misfit = do_mcmc_step(localpart(m)[1], localpart(mn)[1],
        localpart(opt)[1], localpart(stat)[1],localpart(current_misfit)[1], 
        localpart(F)[1], Temp, isample, wp, chain_idx, master_pid)
end

function do_mcmc_step(mn::ModelNuisance, m::Model,
    optn::OptionsNuisance, statn::Stats, current_misfit::Array{Float64,1},
    F::Operator, Temp::Float64, isample::Int, wpn::Writepointers_nuisance, chain_idx::Int, master_pid::Int)
    # this only gets called for moves on the nuisance chain
    # can be either nonstat or stat model
    movetype = do_move!(mn, optn, statn)
    mh_step!(mn, m, F, optn, statn, Temp, movetype, current_misfit)
    get_acceptance_stats!(isample, optn, statn)
    writemodel = false
    abs(Temp - 1.0) < 1e-12 && (writemodel = true)
    write_history(isample, optn, mn, current_misfit[1], statn, wpn, Temp, writemodel, chain_idx, master_pid)
    return current_misfit[1]
end

function do_mcmc_step(mn::DArray{ModelNuisance}, m::DArray{S}, 
    optn::DArray{OptionsNuisance}, statn::DArray{Stats},
    current_misfit::DArray{Array{Float64, 1}},
    F::DArray{x}, Temp::Float64, isample::Int,
    wpn::Writepointers_nuisance, chain_idx::Int, master_pid::Int) where {x<:Operator, S<:Model}
    misfit = do_mcmc_step(localpart(mn)[1], localpart(m)[1], 
            localpart(optn)[1], localpart(statn)[1], localpart(current_misfit)[1],
            localpart(F)[1], localpart(Temp)[1], localpart(isample)[1], wpn, chain_idx, master_pid)
end

function do_mcmc_step(m::ModelStat, mns::ModelNonstat, 
    opt::OptionsStat, optns::OptionsNonstat, stat::Stats,
    current_misfit::Array{Float64, 1}, F::Operator,
    Temp::Float64, isample::Int, wp::Writepointers, chain_idx::Int, master_pid::Int)
    # Stationary GP changes which update nonstationary GP
    movetype, priorviolate = do_move!(m, opt, stat, mns, optns)
    if !priorviolate
        mh_step!(m, mns, F, opt, optns, stat, Temp, movetype, current_misfit)
    end
    get_acceptance_stats!(isample, opt, stat)
    writemodel = false
    abs(Temp-1.0) < 1e-12 && (writemodel = true)
    write_history(isample, opt, m, current_misfit[1], stat, wp, Temp, writemodel, chain_idx, master_pid)
    return current_misfit[1]
end

function do_mcmc_step(m::DArray{ModelStat}, mns::DArray{ModelNonstat}, opt::DArray{OptionsStat}, 
    optns::DArray{OptionsNonstat}, stat::DArray{Stats}, 
    current_misfit::DArray{Array{Float64, 1}}, F::DArray{x}, Temp::Float64, 
    isample::Int, wp::Writepointers, chain_idx::Int, master_pid::Int) where x<:Operator
    misfit = do_mcmc_step(localpart(m)[1], localpart(mns)[1],
        localpart(opt)[1], localpart(optns)[1], localpart(stat)[1], 
        localpart(current_misfit)[1], localpart(F)[1],
        Temp, isample, wp, chain_idx, master_pid)
end

function do_mcmc_step(mns::ModelNonstat, m::ModelStat, 
    optns::OptionsNonstat, statns::Stats, current_misfit::Array{Float64, 1},
    F::Operator, Temp::Float64, isample::Int, wpns::Writepointers, chain_idx::Int, master_pid::Int)
    # purely nonstationary GP moves 
    movetype, priorviolate = do_move!(mns, m, optns, statns)
    if !priorviolate
        mh_step!(mns, m, F, optns, statns, Temp, movetype, current_misfit)
    end
    get_acceptance_stats!(isample, optns, statns)
    writemodel = false
    abs(Temp-1.0) < 1e-12 && (writemodel = true)
    write_history(isample, optns, mns, current_misfit[1], statns, wpns, Temp, writemodel, chain_idx, master_pid)
    return current_misfit[1]
end

function do_mcmc_step(mns::DArray{ModelNonstat}, m::DArray{ModelStat},
    optns::DArray{OptionsNonstat}, statns::DArray{Stats},
    current_misfit::DArray{Array{Float64, 1}}, F::DArray{x}, Temp::Float64, 
    isample::Int, wpns::Writepointers, chain_idx::Int, master_pid::Int) where x<:Operator
    misfit = do_mcmc_step(localpart(mns)[1], localpart(m)[1], 
        localpart(optns)[1], localpart(statns)[1],localpart(current_misfit)[1], 
        localpart(F)[1], Temp, isample, wpns, chain_idx, master_pid)
end

function do_mcmc_step(m::ModelStat, opt::OptionsStat, stat::Stats, current_misfit::Array{Float64, 1},
    F::Operator, Temp::Float64, isample::Int, wp::Writepointers, chain_idx::Int, master_pid::Int)
    # purely stationary GP moves
    movetype, priorviolate = do_move!(m, opt, stat)
    if !priorviolate
        mh_step!(m, F, opt, stat, Temp, movetype, current_misfit)
    end
    get_acceptance_stats!(isample, opt, stat)
    writemodel = false
    abs(Temp-1.0) < 1e-12 && (writemodel = true)
    write_history(isample, opt, m, current_misfit[1], stat, wp, Temp, writemodel, chain_idx, master_pid)
    return current_misfit[1]
end

function do_mcmc_step(m::DArray{ModelStat}, opt::DArray{OptionsStat}, stat::DArray{Stats},
    current_misfit::DArray{Array{Float64, 1}}, F::DArray{x}, Temp::Float64, 
    isample::Int, wp::Writepointers, chain_idx::Int, master_pid::Int) where x<:Operator
    misfit = do_mcmc_step(localpart(m)[1], localpart(opt)[1], 
        localpart(stat)[1],localpart(current_misfit)[1], 
        localpart(F)[1], Temp, isample, wp, chain_idx, master_pid)
end

# function close_history(wp::DArray)
#     for (idx, pid) in enumerate(procs(wp))
#         @sync @spawnat pid close_history(wp[idx])
#     end
# end

function getlastiter(opt_in::Options, chains, idx)
    iterlast = 0
    if opt_in.history_mode=="a"
        # only the value of iterlast at last chain is used
        # as we want to avoid history reads
        if idx == length(chains)
            iterlast = history(opt_in, stat=:iter, chain_idx=idx)[end]
        end
        chains[idx].T = history(opt_in, stat=:T, chain_idx=idx)[end]
    end
    iterlast 
end

function init_out_filenames(opt_in::OptionsStat)
    opt_in.costs_filename = "misfits_"*opt_in.fdataname*"s.bin"
    opt_in.fstar_filename = "models_"*opt_in.fdataname*"s.bin"
    opt_in.x_ftrain_filename = "points_"*opt_in.fdataname*"s.bin"
end

function init_out_filenames(opt_in::OptionsNonstat)
    opt_in.costs_filename = "misfits_"*opt_in.fdataname*"ns.bin"
    opt_in.fstar_filename = "models_"*opt_in.fdataname*"ns.bin"
    opt_in.x_ftrain_filename = "points_"*opt_in.fdataname*"ns.bin"
end

function init_out_filenames(optn_in::OptionsNuisance)
    optn_in.costs_filename = "misfits_"*optn_in.fdataname*"nuisance.bin"
    optn_in.vals_filename = "values_nuisance_"*optn_in.fdataname*".bin"
end

function init_file_pointers_and_darrays(opt_in::OptionsStat,
                            optns_in::OptionsNonstat,
                            optn_in::OptionsNuisance,
                            F_in::Operator, chains::Array{Chain, 1})
    # for nonstat, stat, and nuisances all together                        
    
    init_out_filenames(opt_in)
    init_out_filenames(optns_in)
    init_out_filenames(optn_in)
    
    m_, mns_, mn_, opt_, optns_, optn_, F_in_, stat_, statns_, statn_,
    current_misfit_  = map(x -> Array{Future, 1}(undef, length(chains)), 1:11)
    
    opt_in.history_mode == "a" && setrestartflag.([opt_in, optns_in, optn_in])

    wp = open_history(opt_in)
    wpns = open_history(optns_in)
    wpn = open_history(optn_in)
    
    iterlast = 0
    @sync for(idx, chain) in enumerate(chains)

        opt_[idx]            = @spawnat chain.pid [opt_in]
        optns_[idx]          = @spawnat chain.pid [optns_in]
        optn_[idx]           = @spawnat chain.pid [optn_in]

        m_[idx]              = @spawnat chain.pid [init(opt_in, idx)]
        mns_[idx]            = @spawnat chain.pid [init(optns_in, fetch(m_[idx])[1], idx)]
        mn_[idx]             = @spawnat chain.pid [init(optn_in, idx)]

        stat_[idx]           = @spawnat chain.pid [Stats()]
        statns_[idx]         = @spawnat chain.pid [Stats()]
        statn_[idx]          = @spawnat chain.pid [Stats(nmoves=optn_in.nnu)]

        F_in_[idx]           = @spawnat chain.pid [F_in]

        current_misfit_[idx] = @spawnat chain.pid [[ get_misfit(fetch(mns_[idx])[1],
                                fetch(mn_[idx])[1],
                                fetch(optns_[idx])[1],
                                fetch(F_in_[idx])[1]) ]]
        
        iterlast = getlastiter(opt_in, chains, idx)                        
    end

    m, mns, mn, opt, optns, optn, stat, statns, statn, F,
    current_misfit = map(x -> DArray(x), (m_, mns_, mn_, opt_, optns_, optn_,
                                    stat_, statns_, statn_, F_in_, current_misfit_))
    return m, mns, mn, opt, optns, optn, stat, statns, statn, F, current_misfit,
            wp, wpns, wpn, iterlast
end

function domcmciters(iterlast, nsamples, chains, mns::DArray{ModelNonstat}, m::DArray{ModelStat}, 
            mn::DArray{ModelNuisance}, optns::DArray{OptionsNonstat}, opt::DArray{OptionsStat},
            optn::DArray{OptionsNuisance}, statns, stat, statn, current_misfit, F, wpns, wp, wpn, nominaltime)
    # for nonstat, stat, and nuisances all together         
    
    t, tlong = map(x->time(), 1:2)
    for isample = iterlast+1:iterlast+nsamples
        # we do need each remotecall to finish before 
        # moving on to the next kind of move
        swap_temps(chains)
        @sync for (chain_idx, chain) in enumerate(chains)
            # purely nonstationary GP moves + nuisance
            @async chain.misfit = remotecall_fetch(do_mcmc_step, chain.pid,
                                            mns, m, mn, optns, statns,
                                            current_misfit, F,
                                            chain.T, isample, wpns, chain_idx, chain.master_pid)
        end
        @sync for (chain_idx, chain) in enumerate(chains)
            # purely nuisance move
            @async chain.misfit = remotecall_fetch(do_mcmc_step, chain.pid,
                                            mn, m, optn, statn,
                                            current_misfit, F,
                                            chain.T, isample, wpn, chain_idx, chain.master_pid)
        end
        @sync for (chain_idx, chain) in enumerate(chains)
            # Stationary GP changes which update nonstationary GP + nuisance
            @async chain.misfit = remotecall_fetch(do_mcmc_step, chain.pid,
                                            m, mns, mn, opt, optns, stat,
                                            current_misfit, F,
                                            chain.T, isample, wp, chain_idx, chain.master_pid)
        end
        t, tlong, doquit = disptime(isample, t, tlong, iterlast, nsamples, nominaltime)
        doquit && break
    end
end

function main(opt_in     ::OptionsStat,
            optns_in     ::OptionsNonstat,
            optn_in      ::OptionsNuisance,
            F_in         ::Operator,
            chainprocs   ::Array{Int, 1};
            nsamples     = 4001,
            nchainsatone = 1,
            Tmax         = 2.5,
            nominaltime  = nothing)
    # for nonstat, stat, and nuisances all together 

    chains = Chain(myid(), chainprocs, nchainsatone=nchainsatone, Tmax=Tmax)

    m, mns, mn, opt, optns, optn, stat, 
    statns, statn, F, current_misfit, wp, wpns, wpn, 
    iterlast = init_file_pointers_and_darrays(opt_in, optns_in, optn_in, F_in, chains)

    domcmciters(iterlast, nsamples, chains, mns, m, mn, optns, opt,
        optn, statns, stat, statn, current_misfit, F, wpns, wp, wpn, nominaltime)

    close_history.([wp, wpns, wpn])
    nothing
end

function main(opt_in ::OptionsStat,
    optns_in     ::OptionsNonstat,
    optn         ::OptionsNuisance,
    F_in         ::Operator;
    nchains      = 4,
    nsamples     = 4001,
    nchainsatone = 1,
    Tmax         = 2.5,
    nominaltime  = nothing)
    # unnecessary but for backwards compat
    # for nonstat and stat together 
    pids = workers()
    @assert length(pids) == nchains
    main(opt_in, optns_in, optn, F_in, workers(); nsamples=nsamples, Tmax=Tmax, nchainsatone=nchainsatone, nominaltime)
end

function init_file_pointers_and_darrays(opt_in::OptionsStat,
                            optns_in::OptionsNonstat,
                            F_in::Operator, chains::Array{Chain, 1})
    # for nonstat and stat together                    

    init_out_filenames(opt_in)
    init_out_filenames(optns_in)

    m_, mns_, opt_, optns_, F_in_, stat_, statns_, 
    current_misfit_ = map(x -> Array{Future, 1}(undef, length(chains)), 1:8)
    
    opt_in.history_mode == "a" && setrestartflag.([opt_in, optns_in])

    wp = open_history(opt_in)
    wpns = open_history(optns_in)
    
    iterlast = 0
    @sync for(idx, chain) in enumerate(chains)

        opt_[idx]            = @spawnat chain.pid [opt_in]
        optns_[idx]          = @spawnat chain.pid [optns_in]

        m_[idx]              = @spawnat chain.pid [init(opt_in, idx)]
        mns_[idx]            = @spawnat chain.pid [init(optns_in, fetch(m_[idx])[1], idx)]

        stat_[idx]           = @spawnat chain.pid [Stats()]
        statns_[idx]         = @spawnat chain.pid [Stats()]

        F_in_[idx]           = @spawnat chain.pid [F_in]

        current_misfit_[idx] = @spawnat chain.pid [[ get_misfit(fetch(mns_[idx])[1],
                                fetch(optns_[idx])[1],
                                fetch(F_in_[idx])[1]) ]]
        
        iterlast = getlastiter(opt_in::Options, chains, idx)                        
    end
    
    m, mns, opt, optns, stat, statns, F,
    current_misfit = map(x -> DArray(x), (m_, mns_, opt_, optns_, 
                                    stat_, statns_, F_in_, current_misfit_,
                                    ))
    return m, mns, opt, optns, stat, statns, F, current_misfit,
            wp, wpns, iterlast
end

function domcmciters(iterlast, nsamples, chains, mns::DArray{ModelNonstat}, m::DArray{ModelStat}, 
            optns::DArray{OptionsNonstat}, opt::DArray{OptionsStat}, 
            statns, stat, current_misfit, F, wpns, wp, nominaltime)
    # for nonstat and stat together        
    
    t, tlong = map(x->time(), 1:2)
    for isample = iterlast+1:iterlast+nsamples
        # we do need each remotecall to finish before 
        # moving on to the next kind of move
        swap_temps(chains)
        @sync for (chain_idx, chain) in enumerate(chains)
            # purely nonstationary GP moves
            @async chain.misfit = remotecall_fetch(do_mcmc_step, chain.pid,
                                            mns, m, optns, statns,
                                            current_misfit, F,
                                            chain.T, isample, wpns, chain_idx, chain.master_pid)
        end
        @sync for (chain_idx, chain) in enumerate(chains)
            # Stationary GP changes which update nonstationary GP
            @async chain.misfit = remotecall_fetch(do_mcmc_step, chain.pid,
                                            m, mns, opt, optns, stat,
                                            current_misfit, F,
                                            chain.T, isample, wp, chain_idx, chain.master_pid)
        end
        t, tlong, doquit = disptime(isample, t, tlong, iterlast, nsamples, nominaltime)
        doquit && break
    end
end

function main(opt_in     ::OptionsStat,
            optns_in     ::OptionsNonstat,
            F_in         ::Operator,
            chainprocs   ::Array{Int, 1};
            nsamples     = 4001,
            nchainsatone = 1,
            Tmax         = 2.5,
            nominaltime  = nothing)
    # for nonstat and stat together 

    chains = Chain(myid(), chainprocs, nchainsatone=nchainsatone, Tmax=Tmax)

    m, mns, opt, optns, stat, 
    statns, F, current_misfit, wp, wpns,  
    iterlast = init_file_pointers_and_darrays(opt_in, optns_in, F_in, chains)

    domcmciters(iterlast, nsamples, chains, mns, m, optns, opt,
        statns, stat, current_misfit, F, wpns, wp, nominaltime)

    close_history.([wp, wpns])
    nothing
end

function main(opt_in ::OptionsStat,
    optns_in     ::OptionsNonstat,
    F_in         ::Operator;
    nchains      = 4,
    nsamples     = 4001,
    nchainsatone = 1,
    Tmax         = 2.5,
    nominaltime  = nothing)
    # unnecessary but for backwards compat
    # for nonstat and stat together 
    pids = workers()
    @assert length(pids) == nchains
    main(opt_in, optns_in, F_in, workers(); nsamples=nsamples, Tmax=Tmax, nchainsatone=nchainsatone, nominaltime)
end   

function init_file_pointers_and_darrays(opt_in::OptionsStat, optn_in::OptionsNuisance,
                            F_in::Operator, chains::Array{Chain, 1})
    # purely stationary GP moves + nuisance                    
    
    init_out_filenames(opt_in)
    init_out_filenames(optn_in)

    m_, mn_, opt_, optn_, F_in_, stat_, statn_,
    current_misfit_ = map(x -> Array{Future, 1}(undef, length(chains)), 1:8)
    
    opt_in.history_mode == "a" && setrestartflag.([opt_in, optn_in])

    wp = open_history(opt_in)
    wpn = open_history(optn_in)
    
    iterlast = 0
    @sync for(idx, chain) in enumerate(chains)

        opt_[idx]            = @spawnat chain.pid [opt_in]
        optn_[idx]           = @spawnat chain.pid [optn_in]

        m_[idx]              = @spawnat chain.pid [init(opt_in, idx)]
        mn_[idx]             = @spawnat chain.pid [init(optn_in, idx)]

        stat_[idx]           = @spawnat chain.pid [Stats()]
        statn_[idx]          = @spawnat chain.pid [Stats(nmoves=optn_in.nnu)]

        F_in_[idx]           = @spawnat chain.pid [F_in]    

        current_misfit_[idx] = @spawnat chain.pid [[ get_misfit(fetch(m_[idx])[1],
                                fetch(mn_[idx])[1],
                                fetch(opt_[idx])[1],
                                fetch(F_in_[idx])[1]) ]]
        
        iterlast = getlastiter(opt_in, chains, idx)                         
    end

    m, mn, opt, optn, stat, statn, F,
    current_misfit = map(x -> DArray(x), (m_, mn_, opt_, optn_,
                                    stat_, statn_, F_in_, current_misfit_,
                                    ))
    return m, mn, opt, optn, stat, statn, F, current_misfit,
        wp, wpn, iterlast
end

function domcmciters(iterlast, nsamples, chains, m::DArray{ModelStat}, mn::DArray{ModelNuisance}, 
            opt::DArray{OptionsStat}, optn::DArray{OptionsNuisance}, stat, statn, 
            current_misfit, F, wp, wpn, nominaltime)
    # purely stationary GP moves + nuisance        
    
    t, tlong = map(x->time(), 1:2)
    for isample = iterlast+1:iterlast+nsamples
        # we do need each remotecall to finish before 
        # moving on to the next kind of move
        swap_temps(chains)
        @sync for (chain_idx, chain) in enumerate(chains)
            # purely nuisance move
            @async chain.misfit = remotecall_fetch(do_mcmc_step, chain.pid,
                                            mn, m, optn, statn,
                                            current_misfit, F,
                                            chain.T, isample, wpn, chain_idx, chain.master_pid)
        end
        @sync for (chain_idx, chain) in enumerate(chains)
            # purely stationary GP moves + nuisance
            @async chain.misfit = remotecall_fetch(do_mcmc_step, chain.pid,
                                            m, mn, opt, stat,
                                            current_misfit, F,
                                            chain.T, isample, wp, chain_idx, chain.master_pid)
        end
        t, tlong, doquit = disptime(isample, t, tlong, iterlast, nsamples, nominaltime)
        doquit && break
    end
end

function main(opt_in       ::OptionsStat,
        optn_in      ::OptionsNuisance,
        F_in         ::Operator,
        chainprocs   ::Array{Int, 1};
        nsamples     = 4001,
        nchainsatone = 1,
        Tmax         = 2.5,
        nominaltime  = nothing)
    # purely stationary GP moves + nuisance   

    chains = Chain(myid(), chainprocs, nchainsatone=nchainsatone, Tmax=Tmax)

    m, mn, opt, optn, stat, 
    statn, F, current_misfit, wp, wpn, 
    iterlast = init_file_pointers_and_darrays(opt_in, optn_in, F_in, chains)

    domcmciters(iterlast, nsamples, chains, m, mn, opt,
        optn, stat, statn, current_misfit, F, wp, wpn, nominaltime)

    close_history.([wp, wpn])
    nothing
end

function main(opt_in ::OptionsStat,
    optn_in     ::OptionsNuisance,
    F_in         ::Operator;
    nchains      = 4,
    nsamples     = 4001,
    nchainsatone = 1,
    Tmax         = 2.5,
    nominaltime  = nothing)
    # unnecessary but for backwards compat
    # purely stationary GP moves + nuisance  
    pids = workers()
    @assert length(pids) == nchains
    main(opt_in, optn_in, F_in, workers(); nsamples=nsamples, Tmax=Tmax, nchainsatone=nchainsatone, nominaltime)
end 

function init_file_pointers_and_darrays(opt_in::OptionsStat, F_in::Operator, 
    chains::Array{Chain, 1})
    # purely stationary GP moves
    
    init_out_filenames(opt_in)

    m_, opt_, F_in_, stat_,
    current_misfit_ = map(x -> Array{Future, 1}(undef, length(chains)), 1:5)
    
    opt_in.history_mode == "a" && setrestartflag(opt_in)

    wp = open_history(opt_in)
    
    iterlast = 0
    @sync for(idx, chain) in enumerate(chains)

        opt_[idx]            = @spawnat chain.pid [opt_in]
        m_[idx]              = @spawnat chain.pid [init(opt_in, idx)]
        stat_[idx]           = @spawnat chain.pid [Stats()]
        F_in_[idx]           = @spawnat chain.pid [F_in]    
        current_misfit_[idx] = @spawnat chain.pid [[ get_misfit(fetch(m_[idx])[1],
                                fetch(opt_[idx])[1],
                                fetch(F_in_[idx])[1]) ]]
        iterlast = getlastiter(opt_in, chains, idx)                         
    end

    m, opt, stat, F,
    current_misfit = map(x -> DArray(x), (m_, opt_, 
                                    stat_, F_in_, current_misfit_,
                                    ))
    return m, opt, stat, F, current_misfit,
        wp, iterlast
end

function domcmciters(iterlast, nsamples, chains, m::DArray{ModelStat}, 
            opt::DArray{OptionsStat}, stat, 
            current_misfit, F, wp, nominaltime)
    # purely stationary GP moves     
    
    t, tlong = map(x->time(), 1:2)
    for isample = iterlast+1:iterlast+nsamples
        swap_temps(chains)
        @sync for (chain_idx, chain) in enumerate(chains)
            # purely stationary GP moves 
            @async chain.misfit = remotecall_fetch(do_mcmc_step, chain.pid,
                                            m, opt, stat,
                                            current_misfit, F,
                                            chain.T, isample, wp, chain_idx, chain.master_pid)
        end
        t, tlong, doquit = disptime(isample, t, tlong, iterlast, nsamples, nominaltime)
        doquit && break
    end
end

function main(opt_in ::OptionsStat,
        F_in         ::Operator,
        chainprocs   ::Array{Int, 1};
        nsamples     = 4001,
        nchainsatone = 1,
        Tmax         = 2.5,
        nominaltime  = nothing)
    # purely stationary GP moves    

    chains = Chain(myid(), chainprocs, nchainsatone=nchainsatone, Tmax=Tmax)

    m, opt, stat, 
    F, current_misfit, wp, 
    iterlast = init_file_pointers_and_darrays(opt_in, F_in, chains)

    domcmciters(iterlast, nsamples, chains, m, opt,
        stat, current_misfit, F, wp, nominaltime)

    close_history(wp)
    nothing
end

function main(opt_in ::OptionsStat,
    F_in         ::Operator;
    nchains      = 4,
    nsamples     = 4001,
    nchainsatone = 1,
    Tmax         = 2.5,
    nominaltime  = nothing)
    # unnecessary but for backwards compat
    # purely stationary GP moves    
    pids = workers()
    @assert length(pids) == nchains
    main(opt_in, F_in, workers(); nsamples=nsamples, Tmax=Tmax, nchainsatone=nchainsatone, nominaltime)
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

function disptime(isample, t, tlong, iterlast, nsamples, nominaltime)
    doquit = false
    windowtime = 1000
    nw = 5
    if mod(isample-1, windowtime) == 0
        dt = time() - t #seconds
        t = time()
        @info("on pid $(myid()) **$(@sprintf("%.2f", dt))**sec** $isample out of $(iterlast+nsamples)")
    end
    if !isnothing(nominaltime)
        if mod(isample-1, nw*windowtime) == 0
            dt = time() - tlong #seconds
            tlong = time()
            if (dt/nw) > nominaltime
                doquit = true
                @info("QUIT: PID: $(myid()) TIME:**$(@sprintf("%.2f", dt/nw))**SEC AVG** SAMPLE: $isample")
            end
        end
    end    
    t, tlong, doquit
end

function write_to_log(fp::IOStream, msg)
    write(fp, msg)
    flush(fp)
end

function write_to_bin(fp::IOStream, cidx::Int, data)
    write(fp, cidx)
    write(fp, data)
    flush(fp)
end
