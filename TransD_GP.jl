module TransD_GP
using Printf, LinearAlgebra, Statistics, Distributed, Distances

mutable struct Options
    nmin                :: Int
    nmax                :: Int
    xbounds             :: Array{Float64}
    fbounds             :: Array{Float64}
    xall                :: Array{Float64}
    λ²                 :: Array{Float64,1}
    δ                   :: Float64
    demean              :: Bool
    sdev_prop           :: Float64
    sdev_pos            :: Array{Float64,1}
    pnorm               :: Float64
    stat_window         :: Int
    dispstatstoscreen   :: Bool
    report_freq         :: Int
    save_freq           :: Int
    fdataname           :: String
    history_mode        :: String
    costs_filename      :: String
    fstar_filename      :: String
    x_ftrain_filename   :: String
    debug               :: Bool
    quasimultid         :: Bool
    influenceradius     :: Array{Float64, 1}
end

function Options(;
        nmin               = 1,
        nmax               = 50,
        xbounds            = [0 1.794; 0 1.2],
        fbounds            = [-1 1],
        xall               = nothing,
        λ                 = [0.3],
        δ                 = 0.1,
        demean             = true,
        sdev_prop          = 0.01,
        sdev_pos           = [0.05;0.05],
        pnorm              = 2,
        stat_window        = 100,
        dispstatstoscreen  = true,
        report_freq        = 10,
        save_freq          = 100,
        history_mode       = "w",
        fdataname          = "",
        debug              = false,
        quasimultid        = "",
        influenceradius    = [-9.9]
        )

        @assert xall != nothing
        @assert all(diff(xbounds, dims=2) .> 0)
        @assert diff(fbounds, dims=2)[1] > 0
        @assert ndims(sdev_pos) == 1
        @assert length(sdev_pos) == size(xbounds, 1)
        @assert length(λ) == size(xbounds, 1)
        @assert quasimultid != "" "specify true or false explicitly"
        if quasimultid
            @assert influenceradius[1] > 0.0
            @assert length(influenceradius) == size(xall, 1) - 1
        end
        costs_filename = "misfits_"*fdataname*".bin"
        fstar_filename = "models_"*fdataname*".bin"
        x_ftrain_filename = "points_"*fdataname*".bin"

        Options(nmin, nmax, xbounds, fbounds, xall, λ.^2 , δ, demean, sdev_prop, sdev_pos, pnorm,
                stat_window, dispstatstoscreen, report_freq, save_freq,
                fdataname, history_mode, costs_filename, fstar_filename, x_ftrain_filename,
                debug, quasimultid, influenceradius)
end

mutable struct Model
    fstar         :: Array{Float64}
    xtrain        :: Array{Float64}
    ftrain        :: Array{Float64}
    K_y           :: Array{Float64, 2}
    Kstar         :: Array{Float64, 2}
    n             :: Int
    ftrain_old    :: Float64
    xtrain_old    :: Array{Float64}
    iremember     :: Int # stores old changed point index to recover state
    xtrain_focus  :: Array{Float64} # for quasimultid to keep track of
end

mutable struct Stats
    move_tries::Array{Int, 1}
    accepted_moves::Array{Int, 1}
    accept_rate::Array{Float64, 1}
end

function Stats(;nmoves=4)
    Stats(zeros(Int, nmoves), zeros(Int, nmoves), zeros(Float64, nmoves))
end

mutable struct Writepointers
    fp_costs      :: IOStream
    fp_fstar      :: IOStream
    fp_x_ftrain   :: IOStream
end

function gaussiankernel(x::AbstractArray, y::AbstractArray, p)
    return exp(-0.5*(norm(x-y,p)^p))
end

function init(opt::TransD_GP.Options)
    n = opt.nmin
    xtrain = zeros(Float64, size(opt.xbounds,1), opt.nmax)
    xtrain[:,1:n] = opt.xbounds[:,1] .+ diff(opt.xbounds, dims=2).*rand(size(opt.xbounds, 1), n)
    ftrain = zeros(Float64, opt.nmax)
    ftrain[1:n] = opt.fbounds[1] .+ diff(opt.fbounds, dims=2)[1]*rand(n)
    K_y = zeros(opt.nmax, opt.nmax)
    map!(x²->exp(-x²/2),K_y,pairwise(WeightedSqEuclidean(1 ./opt.λ² ), xtrain, dims=2))
    K_y[diagind(K_y)] .+= opt.δ^2
    Kstar = zeros(Float64, size(opt.xall,2), opt.nmax)
    xtest = opt.xall
    map!(x²->exp(-x²/2),Kstar,pairwise(WeightedSqEuclidean(1 ./opt.λ² ), xtest, xtrain, dims=2))
    mf = 0.
    if opt.demean && n>1
        mf = mean(ftrain[1:n])
    end
    rhs = ftrain[1:n] .- mf
    U = cholesky(K_y[1:n,1:n]).U
    fstar = mf .+ Kstar[:,1:n]*(U\(U'\rhs))
    return Model(fstar, xtrain, ftrain, K_y, Kstar, n,
                 0.0, zeros(Float64, size(opt.xbounds, 1)), 0, zeros(Float64, size(opt.xbounds, 1)))
end

function birth!(m::Model, opt::TransD_GP.Options)
    xtrain, ftrain, K_y,  Kstar, n = m.xtrain, m.ftrain, m.K_y,  m.Kstar, m.n
    xtrain[:,n+1] = opt.xbounds[:,1] + diff(opt.xbounds, dims=2).*rand(size(opt.xbounds, 1))
    copy!(m.xtrain_focus, xtrain[:,n+1])
    ftrain[n+1] = opt.fbounds[1] + diff(opt.fbounds, dims=2)[1]*rand()
    xtest = opt.xall
    Kstarv = @view Kstar[:,n+1]
    map!(x²->exp(-x²/2),Kstarv,colwise(WeightedSqEuclidean(1 ./opt.λ² ), xtrain[:,n+1], xtest))
    K_yv = @view K_y[n+1,1:n+1]
    map!(x²->exp(-x²/2),K_yv,colwise(WeightedSqEuclidean(1 ./opt.λ² ), xtrain[:,n+1], xtrain[:,1:n+1]))
    K_y[1:n+1,n+1] = K_y[n+1,1:n+1]
    K_y[n+1,n+1] = K_y[n+1,n+1] + opt.δ^2
    mf = 0.
    if opt.demean
        mf = mean(ftrain[1:n+1])
    end
    rhs = ftrain[1:n+1] .- mf
    U = cholesky(K_y[1:n+1, 1:n+1]).U
    m.fstar[:] = mf .+ Kstar[:,1:n+1]*(U\(U'\rhs))
    m.n = n+1
end

function undo_birth!(m::Model, opt::TransD_GP.Options)
    m.n = m.n - 1
end

function death!(m::Model, opt::TransD_GP.Options)
    xtrain, ftrain, K_y,  Kstar, n = m.xtrain, m.ftrain, m.K_y,  m.Kstar, m.n
    ipoint = rand(1:n)
    copy!(m.xtrain_focus, xtrain[:,ipoint])
    xtrain[:,ipoint], xtrain[:,n] = xtrain[:,n], xtrain[:,ipoint]
    ftrain[ipoint], ftrain[n] = ftrain[n], ftrain[ipoint]
    Kstar[:,ipoint], Kstar[:,n] = Kstar[:,n], Kstar[:,ipoint]
    K_y[ipoint,1:n], K_y[n,1:n] = K_y[n,1:n], K_y[ipoint,1:n]
    K_y[1:n-1,ipoint] = K_y[ipoint,1:n-1]
    K_y[ipoint,ipoint] = 1.0 + opt.δ^2
    mf = 0.
    if opt.demean && n>2
        mf = mean(ftrain[1:n-1])
    end
    rhs = ftrain[1:n-1] .- mf
    U = cholesky(K_y[1:n-1, 1:n-1]).U
    m.fstar[:] = mf .+ Kstar[:,1:n-1]*(U\(U'\rhs))
    m.n = n-1
end

function undo_death!(m::Model, opt::TransD_GP.Options)
    m.n = m.n + 1
    xtrain, K_y, n = m.xtrain, m.K_y, m.n
    K_yv = @view K_y[n,1:n]
    map!(x²->exp(-x²/2),K_yv,colwise(WeightedSqEuclidean(1 ./opt.λ² ), xtrain[:,n], xtrain[:,1:n]))
    K_y[1:n,n] = K_y[n,1:n]
    K_y[n,n] = K_y[n,n] + opt.δ^2
end

function property_change!(m::Model, opt::TransD_GP.Options)
    ftrain, K_y, Kstar, n = m.ftrain, m.K_y, m. Kstar, m.n
    ipoint = 1 + floor(Int, rand()*n)
    m.iremember = ipoint
    m.ftrain_old = ftrain[ipoint]
    copy!(m.xtrain_focus, m.xtrain[:,ipoint])
    ftrain[ipoint] = ftrain[ipoint] + opt.sdev_prop*randn()
    while (ftrain[ipoint]<opt.fbounds[1]) || (ftrain[ipoint]>opt.fbounds[2])
            (ftrain[ipoint]<opt.fbounds[1]) && (ftrain[ipoint] = 2*opt.fbounds[1] - ftrain[ipoint])
            (ftrain[ipoint]>opt.fbounds[2]) && (ftrain[ipoint] = 2*opt.fbounds[2] - ftrain[ipoint])
    end
    mf = 0.
    if opt.demean
        mf = mean(ftrain[1:n])
    end
    rhs = ftrain[1:n] .- mf
    # could potentially store chol if very time consuming
    U = cholesky(K_y[1:n, 1:n]).U
    m.fstar[:] = mf .+ Kstar[:,1:n]*(U\(U'\rhs))
end

function undo_property_change!(m::Model, opt::TransD_GP.Options)
    ipoint, ftrain = m.iremember, m.ftrain
    ftrain[ipoint] = m.ftrain_old
end

function position_change!(m::Model, opt::TransD_GP.Options)
    xtrain, ftrain, K_y, Kstar, n = m.xtrain, m.ftrain, m.K_y, m.Kstar, m.n
    ipoint = 1 + floor(Int, rand()*n)
    m.iremember = ipoint
    m.xtrain_old[:] = xtrain[:,ipoint]
    xtrain[:,ipoint] = xtrain[:,ipoint] + opt.sdev_pos.*randn(size(opt.xbounds, 1))
    copy!(m.xtrain_focus, xtrain[:,ipoint])
    for i in eachindex(xtrain[:,ipoint])
        while (xtrain[i,ipoint]<opt.xbounds[i,1]) || (xtrain[i,ipoint]>opt.xbounds[i,2])
                (xtrain[i,ipoint]<opt.xbounds[i,1]) && (xtrain[i,ipoint] = 2*opt.xbounds[i,1] - xtrain[i,ipoint])
                (xtrain[i,ipoint]>opt.xbounds[i,2]) && (xtrain[i,ipoint] = 2*opt.xbounds[i,2] - xtrain[i,ipoint])
        end
    end
    xtest = opt.xall
    Kstarv = @view Kstar[:,ipoint]
    map!(x²->exp(-x²/2),Kstarv,colwise(WeightedSqEuclidean(1 ./opt.λ² ), xtrain[:,ipoint], xtest))
    K_yv = @view K_y[ipoint,1:n]
    map!(x²->exp(-x²/2),K_yv,colwise(WeightedSqEuclidean(1 ./opt.λ² ), xtrain[:,ipoint], xtrain[:,1:n]))
    K_y[1:n,ipoint] = K_y[ipoint,1:n]
    K_y[ipoint,ipoint] = K_y[ipoint,ipoint] + opt.δ^2
    mf = 0.
    if opt.demean
        mf = mean(ftrain[1:n])
    end
    rhs = ftrain[1:n] .- mf
    # could potentially store chol if very time consuming
    U = cholesky(K_y[1:n, 1:n]).U
    m.fstar[:] = mf .+ Kstar[:,1:n]*(U\(U'\rhs))
end

function undo_position_change!(m::Model, opt::TransD_GP.Options)
    xtrain, K_y, Kstar, n = m.xtrain, m.K_y, m.Kstar, m.n
    ipoint = m.iremember
    xtrain[:,ipoint] = m.xtrain_old
    xtest = opt.xall
    Kstarv = @view Kstar[:,ipoint]
    map!(x²->exp(-x²/2),Kstarv,colwise(WeightedSqEuclidean(1 ./opt.λ² ), xtrain[:,ipoint], xtest))
    K_yv = @view K_y[ipoint,1:n]
    map!(x²->exp(-x²/2),K_yv,colwise(WeightedSqEuclidean(1 ./opt.λ² ), xtrain[:,ipoint], xtrain[:,1:n]))
    K_y[1:n,ipoint] = K_y[ipoint,1:n]
    K_y[ipoint,ipoint] = K_y[ipoint,ipoint] + opt.δ^2
end

function sync_model!(m::Model, opt::Options)
    ftrain, K_y, Kstar, n = m.ftrain, m.K_y, m.Kstar, m.n
    mf = 0.
    if opt.demean
        mf = mean(ftrain[1:n])
    end
    rhs = ftrain[1:n] .- mf
    # could potentially store chol if very time consuming
    U = cholesky(K_y[1:n, 1:n]).U
    m.fstar[:] = mf .+ Kstar[:,1:n]*(U\(U'\rhs))
end

function do_move!(m::Model, opt::Options, stat::Stats)
    unifrand = rand()
    movetype, priorviolate = 0, false
    if unifrand<0.25
        if m.n<opt.nmax
            birth!(m, opt)
        else
            priorviolate = true
        end
        movetype = 1
    elseif unifrand<0.5
        if m.n>opt.nmin
            death!(m, opt)
        else
            priorviolate = true
        end
        movetype = 2
    elseif unifrand<0.75
        position_change!(m, opt)
        movetype = 3
    else
        property_change!(m, opt)
        movetype = 4
    end
    stat.move_tries[movetype] += 1
    return movetype, priorviolate
end

function undo_move!(movetype::Int, m::Model, opt::Options)
    if movetype == 1
        undo_birth!(m, opt)
    elseif movetype == 2
        undo_death!(m, opt)
    elseif movetype == 3
        undo_position_change!(m, opt)
    else
        undo_property_change!(m, opt)
    end
    sync_model!(m, opt)
    nothing
end

function get_acceptance_stats!(isample::Int, opt::Options, stat::Stats)
    if mod(isample-1, opt.stat_window) == 0
        stat.accept_rate[:] = 100. *stat.accepted_moves./stat.move_tries
        if opt.dispstatstoscreen
            msg = @sprintf("Acceptance rates Birth %5.2f Death %5.2f Position %5.2f Property %5.2f",
                            stat.accept_rate[1],
                            stat.accept_rate[2],
                            stat.accept_rate[3],
                            stat.accept_rate[4])
            @info(msg)
        end
        fill!(stat.move_tries, 0)
        fill!(stat.accepted_moves, 0)
    end
end

# history methods
function open_history(opt::Options)
    if opt.report_freq > 0
        @info("running transD_sampler...")
    end
    fp_costs = nothing
    if length(opt.costs_filename) > 0
        fp_costs = open(opt.costs_filename, opt.history_mode)
    end
    fp_models = nothing
    if length(opt.fstar_filename) > 0
        fp_fstar = open(opt.fstar_filename, opt.history_mode)
    end
    fp_x_ftrain = nothing
    if length(opt.x_ftrain_filename) > 0
        fp_x_ftrain = open(opt.x_ftrain_filename, opt.history_mode)
    end
    return Writepointers(fp_costs, fp_fstar, fp_x_ftrain)
end

function close_history(wp::Writepointers)
    if wp.fp_costs != nothing
        close(wp.fp_costs)
    end
    if wp.fp_fstar != nothing
        close(wp.fp_fstar)
    end
    if wp.fp_x_ftrain != nothing
        close(wp.fp_x_ftrain)
    end
    @info "closed files"
end

function clear_history(opt::Options)
    if length(opt.costs_filename) != 0 && isfile(opt.costs_filename) == true
        rm(opt.costs_filename)
    end
    if length(opt.fstar_filename) != 0 && isfile(opt.fstar_filename) == true
        rm(opt.models_decompr_filename)
    end
    if length(opt.x_ftrain_filename) != 0 && isfile(opt.x_ftrain_filename) == true
        rm(opt.x_ftrain_filename)
    end
end

function mode_history(opt::Options, mode::String)
    @assert mode == "w" || mode == "a"
    opt.history_mode = mode
end

function write_history(isample::Int, opt::Options, m::Model, misfit::Float64,
                        stat::Stats, wp::Writepointers, T::Float64, writemodel::Bool)
    write_history(opt, m.fstar, [m.xtrain' m.ftrain], misfit, stat.accept_rate[1],
                        stat.accept_rate[2], stat.accept_rate[3], stat.accept_rate[4], m.n,
                       isample, wp.fp_costs, wp.fp_fstar, wp.fp_x_ftrain, T, writemodel)
end


function write_history(opt::Options, fstar::AbstractArray, x_ftrain::AbstractArray, U::Float64, acceptanceRateBirth::Float64,
                    acceptanceRateDeath::Float64, acceptanceRatePosition::Float64, acceptanceRateProperty::Float64, nodes::Int,
                    iter::Int, fp_costs::Union{IOStream, Nothing}, fp_fstar::Union{IOStream, Nothing},
                    fp_x_ftrain::Union{IOStream, Nothing}, T::Float64, writemodel::Bool)
    if (mod(iter-1, opt.save_freq) == 0 || iter == 1)
        if fp_costs != nothing
            msg = @sprintf("%d %e %e %e %e %d %e %e\n", iter, acceptanceRateBirth, acceptanceRateDeath,
                                        acceptanceRatePosition, acceptanceRateProperty, nodes, U, T)
            write(fp_costs, msg)
            flush(fp_costs)
        end
        if writemodel
            if fp_fstar != nothing
                write(fp_fstar, convert(Array{eltype(Float64)},fstar))
                flush(fp_fstar)
            end
            if fp_x_ftrain != nothing
                write(fp_x_ftrain, convert(Array{eltype(Float64)},x_ftrain))
                flush(fp_x_ftrain)
            end
        end
    end
end

function history(opt::Options; stat=:U)
    for (statname, el, idx) in ((:iter,                   Int,      1),
                                (:acceptanceRateBirth,    Float64,  2),
                                (:acceptanceRateDeath,    Float64,  3),
                                (:acceptanceRatePosition, Float64,  4),
                                (:acceptanceRateProperty, Float64,  5),
                                (:nodes,                  Int,      6),
                                (:U,                      Float64,  7),
                                (:T,                      Float64,  8))

        if stat == statname
            if length(opt.costs_filename) == 0
                @warn("history, requested $(statname), but you haven't stored this information.")
                return []
            end
            fp_costs = open(opt.costs_filename)
            mark(fp_costs)
            X = Array{el}(undef, countlines(fp_costs))
            reset(fp_costs)
            i = 1
            for line in readlines(fp_costs)
                X[i] = el == Float64 ? parse(Float64,split(line)[idx]) : parse(Int,split(line)[idx])
                i += 1
            end
            return X
        end
    end
    if stat == :fstar
        if length(opt.fstar_filename) == 0
            @warn("history, requested fstar, but you haven't stored this information.")
            return []
        end
        iters, rem = divrem(filesize(opt.fstar_filename), size(opt.xall,2) * sizeof(Float64))
        @assert rem == 0
        fp_models = open(opt.fstar_filename)
        fstar = Array{Array{Float64, 1}}(undef, iters)
        for i = 1:iters
            fstar[i] = zeros(Float64, size(opt.xall,2))
            read!(fp_models, fstar[i])
        end
        return fstar
    end
    if stat == :x_ftrain
        if length(opt.x_ftrain_filename) == 0
            @warn("history, requested x_ftrain, but you haven't stored this information.")
            return []
        end
        iters, rem = divrem(filesize(opt.x_ftrain_filename), opt.nmax * sizeof(Float64)*(1+size(opt.xbounds, 1)))
        @assert rem == 0
        fp_models = open(opt.x_ftrain_filename)
        x_ftrain = Array{Array{Float64,2},1}(undef, iters)
        for i = 1:iters
            x_ftrain[i] = zeros(Float64, opt.nmax, 1+size(opt.xbounds, 1))
            read!(fp_models, x_ftrain[i])
        end
        return x_ftrain
    end
    @warn("history, requested stat: $(stat) is not recognized.")
    return []
end

function closestmultbelow(num::Real, mult::Real)
   r = rem(num, mult)
   if r < 1e-12
       return num
   else
       return num - r
   end
end

end
