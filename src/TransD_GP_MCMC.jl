using Printf, LinearAlgebra, Statistics, Distributed, PositiveFactorizations,
      Distances, NearestNeighbors, DelimitedFiles

using .GP

abstract type Options end

mutable struct OptionsStat <: Options
    nmin                :: Int
    nmax                :: Int
    xbounds             :: Array{Float64, 2}
    fbounds             :: Array{Float64, 2}
    xall                :: Array{Float64, 2}
    λ²                  :: Array{Float64}
    δ                   :: Float64
    demean              :: Bool
    sdev_prop           :: Array{Float64, 1}
    sdev_pos            :: Array{Float64, 1}
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
    K                   :: GP.Kernel
    kdtree              :: KDTree
    balltree            :: BallTree
    timesλ              :: Float64
    needλ²fromlog       :: Bool
    updatenonstat       :: Bool
    updatenuisances     :: Bool
    peskycholesky       :: Bool
end

function OptionsStat(;
        nmin               = 1,
        nmax               = 50,
        xbounds            = [0 1.794; 0 1.2],
        fbounds            = [-1 1],
        xall               = nothing,
        λ                 = [0.3],
        δ                 = 0.1,
        demean             = true,
        sdev_prop          = [0.01],
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
        influenceradius    = [-9.9],
        K                  = SqEuclidean(),
        timesλ             = 1.0,
        needλ²fromlog      = true,
        updatenonstat      = true,
        updatenuisances    = false,
        peskycholesky      = false,
        )

        @assert xall != nothing
        @assert all(diff(xbounds, dims=2) .> 0)
        @assert all(diff(fbounds, dims=2) .> 0)
        @assert ndims(sdev_prop) == 1
        @assert length(sdev_prop) == size(fbounds, 1)
        @assert ndims(sdev_pos) == 1
        @assert length(sdev_pos) == size(xbounds, 1)
        @assert length(λ) == size(xbounds, 1)
        @assert ndims(λ) == 1
        @assert quasimultid != "" "specify true or false explicitly"
        if quasimultid
            @assert influenceradius[1] > 0.0
            @assert length(influenceradius) == size(xall, 1) - 1
        end
        @assert typeof(K) <: GP.Kernel
        costs_filename = "misfits_"*fdataname*".bin"
        fstar_filename = "models_"*fdataname*".bin"
        x_ftrain_filename = "points_"*fdataname*".bin"
        kdtree = KDTree(xall)
        balltree = BallTree(xall./λ)
        OptionsStat(nmin, nmax, xbounds, fbounds, xall, λ.^2 , δ, demean, sdev_prop, sdev_pos, pnorm,
                stat_window, dispstatstoscreen, report_freq, save_freq,
                fdataname, history_mode, costs_filename, fstar_filename, x_ftrain_filename,
                debug, quasimultid, influenceradius, K, kdtree, balltree, timesλ, needλ²fromlog,
                updatenonstat, updatenuisances, peskycholesky)
end

mutable struct OptionsNonstat <: Options
    nmin                :: Int
    nmax                :: Int
    xbounds             :: Array{Float64, 2}
    fbounds             :: Array{Float64, 2}
    xall                :: Array{Float64, 2}
    δ                   :: Float64
    demean              :: Bool
    sdev_prop           :: Array{Float64, 1}
    sdev_pos            :: Array{Float64, 1}
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
    K                   :: GP.Kernel
    kdtree              :: KDTree
    needλ²fromlog       :: Bool
    updatenonstat       :: Bool
    updatenuisances     :: Bool
    peskycholesky       :: Bool
end

function OptionsNonstat(opt::OptionsStat;
        nmin               = 1,
        nmax               = 50,
        fbounds            = [-1 1],
        δ                  = 0.1,
        demean             = true,
        sdev_prop          = [0.01],
        sdev_pos           = [0.05;0.05],
        pnorm              = 2,
        influenceradius    = [-9.9],
        K                  = SqEuclidean()
        )

        @assert all(diff(fbounds, dims=2) .> 0)
        @assert ndims(sdev_prop) == 1
        @assert length(sdev_prop) == size(fbounds, 1)
        @assert ndims(sdev_pos) == 1
        @assert length(sdev_pos) == size(opt.xbounds, 1)
        @assert typeof(K) <: GP.Kernel
        costs_filename = "misfits_ns_"*opt.fdataname*".bin"
        fstar_filename = "models_ns_"*opt.fdataname*".bin"
        x_ftrain_filename = "points_ns_"*opt.fdataname*".bin"
        OptionsNonstat(nmin, nmax, opt.xbounds, fbounds, opt.xall, δ, demean, sdev_prop, sdev_pos, pnorm,
                opt.stat_window, opt.dispstatstoscreen, opt.report_freq, opt.save_freq,
                opt.fdataname, opt.history_mode, opt.costs_filename, opt.fstar_filename, opt.x_ftrain_filename,
                opt.debug, opt.quasimultid, influenceradius, K, opt.kdtree, opt.needλ²fromlog, opt.updatenonstat,
                opt.updatenuisances, opt.peskycholesky)
end

#not subtyping Options because methods that expect Options are expecting
#options for a Gaussian process model, not fixed-dimension MCMC
mutable struct OptionsNuisance
    sdev                   :: Array{Float64,1}
    bounds                 :: Array{Float64,2}
    nnu                    :: Int64

    updatenuisances         :: Bool
    updatenonstat          :: Bool
    debug                  :: Bool

    stat_window            :: Int
    dispstatstoscreen      :: Bool

    report_freq            :: Int
    save_freq              :: Int
    history_mode           :: String
    costs_filename         :: String
    vals_filename          :: String
end

#make an empty nuisance options struct
OptionsNuisance(opt::OptionsStat, updatenonstat=false) =
    OptionsNuisance([], Array{Float64,2}(undef,0,2), 0, false,
     updatenonstat, false, opt.stat_window, opt.dispstatstoscreen, opt.report_freq, opt.save_freq, opt.history_mode, "misfits_nuisance_"*opt.fdataname*".bin",
     "values_nuisance_"*opt.fdataname*".bin")

#"Model" means a Gaussian-process parametrised function
abstract type Model end

mutable struct ModelStat <:Model
    fstar         :: Array{Float64}
    xtrain        :: Array{Float64}
    ftrain        :: Array{Float64}
    K_y           :: Array{Float64, 2}
    Kstar         :: Array{Float64, 2}
    n             :: Int
    ftrain_old    :: Array{Float64}
    xtrain_old    :: Array{Float64}
    iremember     :: Int # stores old changed point index to recover state
    xtrain_focus  :: Array{Float64} # for quasimultid to keep track of
end

mutable struct ModelNonstat <:Model
    fstar         :: Array{Float64}
    xtrain        :: Array{Float64}
    ftrain        :: Array{Float64}
    K_y           :: Array{Float64, 2}
    Kstar         :: Array{Float64, 2}
    n             :: Int
    ftrain_old    :: Array{Float64}
    xtrain_old    :: Array{Float64}
    iremember     :: Int # stores old changed point index to recover state
    xtrain_focus  :: Array{Float64} # for quasimultid to keep track of
    K_y_old       :: Array{Float64, 2}
    Kstar_old     :: Array{Float64, 2}
end

mutable struct ModelNuisance
    #"nuisances" are a fixed-dimensional part of the model, with the only
    #allowed move being a value change. Contrasting with the GP parametrisation
    #of the main part of the model which allows birth and death of GP nodes
    nuisance      :: Array{Float64}
    nidx_mem      :: Int #last nuisance val changed
    nu_old        :: Float64 #prev val of nidx_mem
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

mutable struct Writepointers_nuisance
    fp_costs      :: IOStream
    fp_vals       :: IOStream
end

# Stationary GP functions, i.e., for λ
function init(opt::OptionsStat)
    n, xtrain, ftrain = initvalues(opt)
    K_y = zeros(opt.nmax, opt.nmax)
    map!(x->GP.κ(opt.K, x),K_y,pairwise(WeightedEuclidean(1 ./opt.λ² ), xtrain, dims=2))
    K_y[diagind(K_y)] .+= opt.δ^2
    Kstar = zeros(Float64, size(opt.xall,2), opt.nmax)
    xtest = opt.xall
    map!(x->GP.κ(opt.K, x),Kstar,pairwise(WeightedEuclidean(1 ./opt.λ² ), xtest, xtrain, dims=2))
    fstar = zeros(size(opt.fbounds, 1), size(opt.xall, 2))
    calcfstar!(fstar, ftrain, opt, K_y, Kstar, n)
    return ModelStat(fstar, xtrain, ftrain, K_y, Kstar, n,
                 [0.0], zeros(Float64, size(opt.xbounds, 1)), 0,
                 zeros(Float64, size(opt.xbounds, 1)))
end

function init(opt::OptionsNuisance)
    @info opt
    if length(opt.bounds) > 0
        nuisance = opt.bounds[:,1] + diff(opt.bounds, dims = 2)[:].*rand(opt.nnu)
    else
        nuisance = []
    end
    #TODO load history if specified (like with GP init)
    return ModelNuisance(nuisance, 0, 0.)
end

function calcfstar!(fstar::Array{Float64,2}, ftrain::Array{Float64,2},
                    opt::OptionsStat, K_y::Array{Float64,2},
                    Kstar::Array{Float64, 2}, n::Int)
    if opt.demean && n>1
        mf = mean(ftrain[:,1:n], dims=2)
    else
        mf = mean(opt.fbounds, dims=2)
    end
    rhs = ftrain[:,1:n] .- mf
    ky = @view K_y[1:n,1:n]
    if opt.peskycholesky
        U = cholesky(Positive, ky, Val{false}).U
    else
        U = cholesky(ky).U
    end
    ks = @view Kstar[:,1:n]
    if opt.needλ²fromlog
        copy!(fstar, 10 .^(2(mf' .+ ks*(U\(U'\rhs'))))' )
    else
        copy!(fstar, (mf' .+ ks*(U\(U'\rhs')))' )
    end
    nothing
end

function initvalues(opt::Options)
    xtrain = zeros(Float64, size(opt.xbounds,1), opt.nmax)
    ftrain = zeros(Float64, size(opt.fbounds,1), opt.nmax)
    if opt.history_mode == "w" # new start
        n = opt.nmin
        xtrain[:,1:n] = opt.xbounds[:,1] .+ diff(opt.xbounds, dims=2).*rand(size(opt.xbounds, 1), n)
        ftrain[:,1:n] = opt.fbounds[:,1] .+ diff(opt.fbounds, dims=2).*rand(size(opt.fbounds, 1), n)
    else # restart
        @info "opening $(opt.x_ftrain_filename)"
        n = history(opt, stat=:nodes)[end]
        xft = history(opt, stat=:x_ftrain)[end]
        xtrain[:,1:n] = xft[1:size(opt.xall, 1), 1:n]
        ftrain[:,1:n] = xft[size(opt.xall, 1)+1:end, 1:n]
    end
    n, xtrain, ftrain
end

function updatenskernels!(opt::OptionsStat, m::ModelStat, ipoint::Union{Int, Array{Int, 1}},
                          optns::OptionsNonstat, mns::ModelNonstat; doall=false, isposchange=false)
    xtrain, λ², xtest = m.xtrain, m.fstar, opt.xall
    copyto!(mns.K_y_old, CartesianIndices((mns.n, mns.n)), mns.K_y, CartesianIndices((mns.n, mns.n)))
    copyto!(mns.Kstar_old, CartesianIndices((size(mns.Kstar_old, 1), mns.n)),
            mns.Kstar, CartesianIndices((size(mns.Kstar_old, 1), mns.n)))
    xt = xtrain[:,ipoint]
    isposchange && (xt = hcat(xt, m.xtrain_old))

    if doall
        kstarchangeidx = 1:size(mns.Kstar, 1)
        kychangeidx = 1:mns.n
    else
        kstarchangeidx = unique(vcat(inrange(opt.balltree, xt./sqrt.(opt.λ²), opt.timesλ)...))
        balltree = BallTree(mns.xtrain[:,1:mns.n]./sqrt.(opt.λ²))
        kychangeidx = inrange(balltree, xt./sqrt.(opt.λ²), opt.timesλ)
    end
    idxs = gettrainidx(opt.kdtree, mns.xtrain, mns.n)
    # changes where test points are in influence radius of lscale change
    ks = view(mns.Kstar, kstarchangeidx, 1:mns.n)
    λ²train = @view λ²[:,idxs]
    λ²test = @view λ²[:,kstarchangeidx]
    xtr = @view mns.xtrain[:,1:mns.n]
    xte = @view xtest[:,kstarchangeidx]
    GP.pairwise(ks, optns.K, xtr, xte, λ²train, λ²test)
    if length(kychangeidx) > 0
        isposchange && (kychangeidx = reduce(vcat, kychangeidx))
        # set 1 of changes where training points are in influence radius of lscale change
        ks2 = view(mns.Kstar, :, kychangeidx)
        xtr2 = @view mns.xtrain[:,kychangeidx]
        λ²train2 = @view λ²[:,idxs][:,kychangeidx]
        GP.pairwise(ks2, optns.K, xtr2, xtest, λ²train2, λ²)
        # set 2 of changes where training points are in influence radius of lscale change
        ky = view(mns.K_y, 1:mns.n , kychangeidx)
        GP.pairwise(ky, optns.K, xtr2, xtr, λ²train2, λ²train)
        mns.K_y[kychangeidx,1:mns.n] .= mns.K_y[1:mns.n,kychangeidx]'
        # nugget add
        ky = view(mns.K_y, 1:mns.n , 1:mns.n)
        ky[diagind(ky)] .= 1. + optns.δ^2
    end
    sync_model!(mns, optns)
    nothing
end

function birth!(m::ModelStat, opt::OptionsStat,
                mns::ModelNonstat, optns::OptionsNonstat; doall=false)
    xtrain, ftrain, K_y,  Kstar, n = m.xtrain, m.ftrain, m.K_y,  m.Kstar, m.n
    xtrain[:,n+1] = opt.xbounds[:,1] + diff(opt.xbounds, dims=2).*rand(size(opt.xbounds, 1))
    copy!(m.xtrain_focus, xtrain[:,n+1])
    ftrain[:,n+1] = opt.fbounds[:,1] + diff(opt.fbounds, dims=2).*rand(size(opt.fbounds, 1))
    xtest = opt.xall
    Kstarv = @view Kstar[:,n+1]
    map!(x->GP.κ(opt.K, x),Kstarv,colwise(WeightedEuclidean(1 ./opt.λ² ), xtrain[:,n+1], xtest))
    K_yv = @view K_y[n+1,1:n+1]
    map!(x->GP.κ(opt.K, x),K_yv,colwise(WeightedEuclidean(1 ./opt.λ² ), xtrain[:,n+1], xtrain[:,1:n+1]))
    K_y[1:n+1,n+1] = K_y[n+1,1:n+1]
    K_y[n+1,n+1] = K_y[n+1,n+1] + opt.δ^2
    calcfstar!(m.fstar, m.ftrain, opt, K_y, Kstar, n+1)
    m.n = n+1
    # updating the nonstationary kernels now
    opt.updatenonstat && updatenskernels!(opt, m, n+1, optns, mns, doall=doall)
    nothing
end

function undo_birth!(m::ModelStat, opt::OptionsStat, mns::ModelNonstat)
    m.n = m.n - 1
    # updating the nonstationary kernels now
    if opt.updatenonstat
        copyto!(mns.K_y, CartesianIndices((mns.n, mns.n)), mns.K_y_old, CartesianIndices((mns.n, mns.n)))
        copyto!(mns.Kstar, CartesianIndices((size(mns.Kstar_old, 1), mns.n)),
                mns.Kstar_old, CartesianIndices((size(mns.Kstar_old, 1), mns.n)))
    end
    nothing
end

function death!(m::ModelStat, opt::OptionsStat,
                mns::ModelNonstat, optns::OptionsNonstat; doall=false)
    xtrain, ftrain, K_y,  Kstar, n = m.xtrain, m.ftrain, m.K_y,  m.Kstar, m.n
    ipoint = rand(1:n)
    m.iremember = ipoint
    copy!(m.xtrain_focus, xtrain[:,ipoint])
    xtrain[:,ipoint], xtrain[:,n] = xtrain[:,n], xtrain[:,ipoint]
    ftrain[:,ipoint], ftrain[:,n] = ftrain[:,n], ftrain[:,ipoint]
    Kstar[:,ipoint], Kstar[:,n] = Kstar[:,n], Kstar[:,ipoint]
    K_y[ipoint,1:n], K_y[n,1:n] = K_y[n,1:n], K_y[ipoint,1:n]
    K_y[1:n-1,ipoint] = K_y[ipoint,1:n-1]
    K_y[ipoint,ipoint] = 1.0 + opt.δ^2
    calcfstar!(m.fstar, m.ftrain, opt, K_y, Kstar, n-1)
    m.n = n-1
    # updating the nonstationary kernels now
    opt.updatenonstat && updatenskernels!(opt, m, n, optns, mns, doall=doall)
    nothing
end

function undo_death!(m::ModelStat, opt::OptionsStat, mns::ModelNonstat)
    m.n = m.n + 1
    ftrain, xtrain, Kstar, K_y, n, ipoint = m.ftrain, m.xtrain, m.Kstar, m.K_y, m.n, m.iremember
    xtrain[:,ipoint], xtrain[:,n] = xtrain[:,n], xtrain[:,ipoint]
    ftrain[:,ipoint], ftrain[:,n] = ftrain[:,n], ftrain[:,ipoint]
    Kstar[:,ipoint], Kstar[:,n] = Kstar[:,n], Kstar[:,ipoint]
    K_y[ipoint,1:n], K_y[n,1:n] = K_y[n,1:n], K_y[ipoint,1:n]
    K_y[1:n,ipoint] = K_y[ipoint,1:n]
    K_y[ipoint,ipoint] = 1.0 + opt.δ^2
    K_y[n,n] = 1.0 + opt.δ^2
    # updating the nonstationary kernels now
    if opt.updatenonstat
        copyto!(mns.K_y, CartesianIndices((mns.n, mns.n)), mns.K_y_old, CartesianIndices((mns.n, mns.n)))
        copyto!(mns.Kstar, CartesianIndices((size(mns.Kstar_old, 1), mns.n)),
                mns.Kstar_old, CartesianIndices((size(mns.Kstar_old, 1), mns.n)))
    end
    nothing
end

function property_change!(m::ModelStat, opt::OptionsStat,
                          mns::ModelNonstat, optns::OptionsNonstat; doall=false)
    ftrain, K_y, Kstar, n = m.ftrain, m.K_y, m. Kstar, m.n
    ipoint = 1 + floor(Int, rand()*n)
    m.iremember = ipoint
    copy!(m.ftrain_old, ftrain[:,ipoint])
    ftrain[:,ipoint] = ftrain[:,ipoint] + opt.sdev_prop.*randn(size(opt.fbounds, 1))
    for i in eachindex(ftrain[:,ipoint])
        while (ftrain[i,ipoint]<opt.fbounds[i,1]) || (ftrain[i,ipoint]>opt.fbounds[i,2])
                (ftrain[i,ipoint]<opt.fbounds[i,1]) && (ftrain[i,ipoint] = 2*opt.fbounds[i,1] - ftrain[i,ipoint])
                (ftrain[i,ipoint]>opt.fbounds[i,2]) && (ftrain[i,ipoint] = 2*opt.fbounds[i,2] - ftrain[i,ipoint])
        end
    end
    calcfstar!(m.fstar, m.ftrain, opt, K_y, Kstar, n)
    # updating the nonstationary kernels now
    opt.updatenonstat && updatenskernels!(opt, m, ipoint, optns, mns, doall=doall)
    nothing
end

function undo_property_change!(m::ModelStat, opt::OptionsStat, mns::ModelNonstat)
    ipoint, ftrain = m.iremember, m.ftrain
    ftrain[:,ipoint] = m.ftrain_old
    # updating the nonstationary kernels now
    if opt.updatenonstat
        copyto!(mns.K_y, CartesianIndices((mns.n, mns.n)), mns.K_y_old, CartesianIndices((mns.n, mns.n)))
        copyto!(mns.Kstar, CartesianIndices((size(mns.Kstar_old, 1), mns.n)),
                mns.Kstar_old, CartesianIndices((size(mns.Kstar_old, 1), mns.n)))
    end
    nothing
end

function position_change!(m::ModelStat, opt::OptionsStat,
                          mns::ModelNonstat, optns::OptionsNonstat; doall=false)
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
    map!(x->GP.κ(opt.K, x),Kstarv,colwise(WeightedEuclidean(1 ./opt.λ² ), xtrain[:,ipoint], xtest))
    K_yv = @view K_y[ipoint,1:n]
    map!(x->GP.κ(opt.K, x),K_yv,colwise(WeightedEuclidean(1 ./opt.λ² ), xtrain[:,ipoint], xtrain[:,1:n]))
    K_y[1:n,ipoint] = K_y[ipoint,1:n]
    K_y[ipoint,ipoint] = K_y[ipoint,ipoint] + opt.δ^2
    calcfstar!(m.fstar, m.ftrain, opt, K_y, Kstar, n)
    # updating the nonstationary kernels now
    opt.updatenonstat && updatenskernels!(opt, m, ipoint, optns, mns, doall=doall, isposchange=true)
    nothing
end

function undo_position_change!(m::ModelStat, opt::OptionsStat, mns::ModelNonstat)
    xtrain, K_y, Kstar, n = m.xtrain, m.K_y, m.Kstar, m.n
    ipoint = m.iremember
    xtrain[:,ipoint] = m.xtrain_old
    xtest = opt.xall
    Kstarv = @view Kstar[:,ipoint]
    map!(x->GP.κ(opt.K, x),Kstarv,colwise(WeightedEuclidean(1 ./opt.λ² ), xtrain[:,ipoint], xtest))
    K_yv = @view K_y[ipoint,1:n]
    map!(x->GP.κ(opt.K, x),K_yv,colwise(WeightedEuclidean(1 ./opt.λ² ), xtrain[:,ipoint], xtrain[:,1:n]))
    K_y[1:n,ipoint] = K_y[ipoint,1:n]
    K_y[ipoint,ipoint] = K_y[ipoint,ipoint] + opt.δ^2
    # updating the nonstationary kernels now
    if opt.updatenonstat
        copyto!(mns.K_y, CartesianIndices((mns.n, mns.n)), mns.K_y_old, CartesianIndices((mns.n, mns.n)))
        copyto!(mns.Kstar, CartesianIndices((size(mns.Kstar_old, 1), mns.n)),
                mns.Kstar_old, CartesianIndices((size(mns.Kstar_old, 1), mns.n)))
    end
    nothing
end

# Non stationary GP functions, i.e., for mns.fstar
function gettrainidx(kdtree::KDTree, xtrain::Array{Float64, 2}, n::Int)
    idxs,  = knn(kdtree, xtrain[:,1:n], 1)
    reduce(vcat, idxs)
end

function init(opt::OptionsNonstat, m::ModelStat)
    donotinit = !opt.needλ²fromlog && !opt.updatenonstat
    if !donotinit
        λ² = m.fstar
        n, xtrain, ftrain = initvalues(opt)
        K_y = zeros(opt.nmax, opt.nmax)
        idxs = gettrainidx(opt.kdtree, xtrain, n)
        ky = view(K_y, 1:n, 1:n)
        GP.pairwise(ky, opt.K, xtrain[:,1:n], xtrain[:,1:n], λ²[:,idxs], λ²[:,idxs])
        K_y[diagind(K_y)] .+= opt.δ^2
        Kstar = zeros(Float64, size(opt.xall,2), opt.nmax)
        xtest = opt.xall
        ks = view(Kstar, :, 1:n)
        GP.pairwise(ks, opt.K, xtrain[:,1:n], xtest, λ²[:,idxs], λ²)
        fstar = zeros(size(opt.xall, 2), size(opt.fbounds, 1))
        calcfstar!(fstar, ftrain, opt, K_y, Kstar, n)
        return ModelNonstat(fstar, xtrain, ftrain, K_y, Kstar, n,
                 [0.0], zeros(Float64, size(opt.xbounds, 1)), 0, zeros(Float64, size(opt.xbounds, 1)),
                 copy(K_y), copy(Kstar))
    else
        dummy2d = [0.0 0.0]
        return ModelNonstat(dummy2d, dummy2d, dummy2d, dummy2d, dummy2d, 0,
                 dummy2d, dummy2d, 0, dummy2d,
                 dummy2d, dummy2d)
    end
end

function calcfstar!(fstar::Array{Float64,2}, ftrain::Array{Float64,2},
                    opt::OptionsNonstat, K_y::Array{Float64,2},
                    Kstar::Array{Float64, 2}, n::Int)
    if opt.demean && n>1
        mf = mean(ftrain[:,1:n], dims=2)
    else
        mf = mean(opt.fbounds, dims=2)
    end
    rhs = ftrain[:,1:n] .- mf
    ky = view(K_y, 1:n, 1:n)
    ks = @view Kstar[:,1:n]
    if opt.peskycholesky
        U = cholesky(Positive, ky, Val{false}).U
    else
        U = cholesky(ky).U
    end
    copy!(fstar, mf' .+ ks*(U\(U'\rhs')) )
    nothing
end

function birth!(m::ModelNonstat, opt::OptionsNonstat, l::ModelStat)
    λ² = l.fstar
    xtrain, ftrain, K_y,  Kstar, n = m.xtrain, m.ftrain, m.K_y,  m.Kstar, m.n
    xtrain[:,n+1] = opt.xbounds[:,1] + diff(opt.xbounds, dims=2).*rand(size(opt.xbounds, 1))
    copy!(m.xtrain_focus, xtrain[:,n+1])
    ftrain[:,n+1] = opt.fbounds[:,1] + diff(opt.fbounds, dims=2).*rand(size(opt.fbounds, 1))
    xtest = opt.xall
    Kstarv = @view Kstar[:,n+1]
    idxs = gettrainidx(opt.kdtree, xtrain, n+1)
    K_yv = @view K_y[n+1,1:n+1]
    @views begin
    GP.colwise!(Kstarv, opt.K, xtrain[:,n+1], xtest, λ²[:,idxs[end]], λ²)
    GP.colwise!(K_yv, opt.K, xtrain[:,n+1], xtrain[:,1:n+1], λ²[:,idxs[end]], λ²[:,idxs])
    end
    K_y[1:n+1,n+1] = K_y[n+1,1:n+1]
    K_y[n+1,n+1] = K_y[n+1,n+1] + opt.δ^2
    calcfstar!(m.fstar, m.ftrain, opt, K_y, Kstar, n+1)
    m.n = n+1
    nothing
end

function undo_birth!(m::ModelNonstat)
    m.n = m.n - 1
    nothing
end

function death!(m::ModelNonstat, opt::OptionsNonstat)
    xtrain, ftrain, K_y,  Kstar, n = m.xtrain, m.ftrain, m.K_y,  m.Kstar, m.n
    ipoint = rand(1:n)
    m.iremember = ipoint
    copy!(m.xtrain_focus, xtrain[:,ipoint])
    xtrain[:,ipoint], xtrain[:,n] = xtrain[:,n], xtrain[:,ipoint]
    ftrain[:,ipoint], ftrain[:,n] = ftrain[:,n], ftrain[:,ipoint]
    Kstar[:,ipoint], Kstar[:,n] = Kstar[:,n], Kstar[:,ipoint]
    K_y[ipoint,1:n], K_y[n,1:n] = K_y[n,1:n], K_y[ipoint,1:n]
    K_y[1:n-1,ipoint] = K_y[ipoint,1:n-1]
    K_y[ipoint,ipoint] = 1.0 + opt.δ^2
    calcfstar!(m.fstar, m.ftrain, opt, K_y, Kstar, n-1)
    m.n = n-1
    nothing
end

function undo_death!(m::ModelNonstat, opt::OptionsNonstat)
    m.n = m.n + 1
    ftrain, xtrain, Kstar, K_y, n, ipoint = m.ftrain, m.xtrain, m.Kstar, m.K_y, m.n, m.iremember
    xtrain[:,ipoint], xtrain[:,n] = xtrain[:,n], xtrain[:,ipoint]
    ftrain[:,ipoint], ftrain[:,n] = ftrain[:,n], ftrain[:,ipoint]
    Kstar[:,ipoint], Kstar[:,n] = Kstar[:,n], Kstar[:,ipoint]
    K_y[ipoint,1:n], K_y[n,1:n] = K_y[n,1:n], K_y[ipoint,1:n]
    K_y[1:n,ipoint] = K_y[ipoint,1:n]
    K_y[ipoint,ipoint] = 1.0 + opt.δ^2
    K_y[n,n] = 1.0 + opt.δ^2
end

function property_change!(m::ModelNonstat, opt::OptionsNonstat)
    ftrain, K_y, Kstar, n = m.ftrain, m.K_y, m. Kstar, m.n
    ipoint = 1 + floor(Int, rand()*n)
    m.iremember = ipoint
    copy!(m.ftrain_old, ftrain[:,ipoint])
    ftrain[:,ipoint] = ftrain[:,ipoint] + opt.sdev_prop.*randn(size(opt.fbounds, 1))
    for i in eachindex(ftrain[:,ipoint])
        while (ftrain[i,ipoint]<opt.fbounds[i,1]) || (ftrain[i,ipoint]>opt.fbounds[i,2])
                (ftrain[i,ipoint]<opt.fbounds[i,1]) && (ftrain[i,ipoint] = 2*opt.fbounds[i,1] - ftrain[i,ipoint])
                (ftrain[i,ipoint]>opt.fbounds[i,2]) && (ftrain[i,ipoint] = 2*opt.fbounds[i,2] - ftrain[i,ipoint])
        end
    end
    calcfstar!(m.fstar, m.ftrain, opt, K_y, Kstar, n)
    nothing
end

function undo_property_change!(m::ModelNonstat)
    ipoint, ftrain = m.iremember, m.ftrain
    ftrain[:,ipoint] = m.ftrain_old
    nothing
end

function position_change!(m::ModelNonstat, opt::OptionsNonstat, l::ModelStat)
    λ² = l.fstar
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
    idxs = gettrainidx(opt.kdtree, xtrain, n)
    K_yv = @view K_y[ipoint,1:n]
    @views begin
    GP.colwise!(Kstarv, opt.K, xtrain[:,ipoint], xtest, λ²[:,idxs[ipoint]], λ²)
    GP.colwise!(K_yv, opt.K, xtrain[:,ipoint], xtrain[:,1:n], λ²[:,idxs[ipoint]], λ²[:,idxs])
    end
    K_y[1:n,ipoint] = K_y[ipoint,1:n]
    K_y[ipoint,ipoint] = K_y[ipoint,ipoint] + opt.δ^2
    calcfstar!(m.fstar, m.ftrain, opt, K_y, Kstar, n)
    nothing
end

function undo_position_change!(m::ModelNonstat, opt::OptionsNonstat, l::ModelStat)
    λ² = l.fstar
    xtrain, K_y, Kstar, n = m.xtrain, m.K_y, m.Kstar, m.n
    ipoint = m.iremember
    xtrain[:,ipoint] = m.xtrain_old
    xtest = opt.xall
    Kstarv = @view Kstar[:,ipoint]
    idxs = gettrainidx(opt.kdtree, xtrain, n)
    K_yv = @view K_y[ipoint,1:n]
    @views begin
    GP.colwise!(Kstarv, opt.K, xtrain[:,ipoint], xtest, λ²[:,idxs[ipoint]], λ²)
    GP.colwise!(K_yv, opt.K, xtrain[:,ipoint], xtrain[:,1:n], λ²[:,idxs[ipoint]], λ²[:,idxs])
    end
    K_y[1:n,ipoint] = K_y[ipoint,1:n]
    K_y[ipoint,ipoint] = K_y[ipoint,ipoint] + opt.δ^2
    nothing
end

function sync_model!(m::Model, opt::Options)
    ftrain, K_y, Kstar, n = m.ftrain, m.K_y, m.Kstar, m.n
    calcfstar!(m.fstar, m.ftrain, opt, K_y, Kstar, n)
    nothing
end

function testupdate(optns::OptionsNonstat, l::ModelStat, mns::ModelNonstat)
    λ² = l.fstar
    idxs = gettrainidx(optns.kdtree, mns.xtrain, mns.n)
    ftest, = GP.GPfit(optns.K, mns.ftrain[:,1:mns.n], mns.xtrain[:,1:mns.n],
            optns.xall, λ², λ²[:,idxs], optns.δ, p=2, demean=optns.demean, nogetvars=true,
            my=mean(optns.fbounds, dims=2))
    return ftest
end

function testupdate(opt::OptionsStat, m::ModelStat)
    ftest, = GP.GPfit(opt.K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n],
            opt.xall, opt.λ², opt.δ, nogetvars=true, demean=opt.demean, p=2,
            my=mean(opt.fbounds, dims=2))
    return ftest
end

function do_move!(mns::ModelNonstat, m::ModelStat, optns::OptionsNonstat, statns::Stats)
    unifrand = rand()
    movetype, priorviolate = 0, false
    if unifrand<0.25
        if mns.n<optns.nmax
            birth!(mns, optns, m)
        else
            priorviolate = true
        end
        movetype = 1
    elseif unifrand<0.5
        if mns.n>optns.nmin
            death!(mns, optns)
        else
            priorviolate = true
        end
        movetype = 2
    elseif unifrand<0.75
        position_change!(mns, optns, m)
        movetype = 3
    else
        property_change!(mns, optns)
        movetype = 4
    end
    statns.move_tries[movetype] += 1
    return movetype, priorviolate
end

function undo_move!(movetype::Int, mns::ModelNonstat, optns::OptionsNonstat, m::ModelStat)
    if movetype == 1
        undo_birth!(mns)
    elseif movetype == 2
        undo_death!(mns, optns)
    elseif movetype == 3
        undo_position_change!(mns, optns, m)
    else
        undo_property_change!(mns)
    end
    sync_model!(mns, optns)
    nothing
end

function do_move!(m::ModelStat, opt::OptionsStat, stat::Stats,
                  mns::ModelNonstat, optns::OptionsNonstat)
    unifrand = rand()
    movetype, priorviolate = 0, false
    if unifrand<0.25
        if m.n<opt.nmax
            birth!(m, opt, mns, optns)
        else
            priorviolate = true
        end
        movetype = 1
    elseif unifrand<0.5
        if m.n>opt.nmin
            death!(m, opt, mns, optns)
        else
            priorviolate = true
        end
        movetype = 2
    elseif unifrand<0.75
        position_change!(m, opt, mns, optns)
        movetype = 3
    else
        property_change!(m, opt, mns, optns)
        movetype = 4
    end
    stat.move_tries[movetype] += 1
    return movetype, priorviolate
end

function undo_move!(movetype::Int, m::ModelStat, opt::OptionsStat,
                    mns::ModelNonstat, optns::OptionsNonstat)
    if movetype == 1
        undo_birth!(m, opt, mns)
    elseif movetype == 2
        undo_death!(m, opt, mns)
    elseif movetype == 3
        undo_position_change!(m, opt, mns)
    else
        undo_property_change!(m, opt, mns)
    end
    sync_model!(m, opt)
    opt.updatenonstat && sync_model!(mns, optns)
    nothing
end

function do_move!(mn::ModelNuisance, optn::OptionsNuisance, statn::Stats)
    priorviolate = false
    nidx = 0
    while true
        nidx = rand(1:optn.nnu)
        numin, numax = extrema(optn.bounds[nidx,:])
        if abs(numin-numax)>1e-12
             break
        end
    end
    new_nval = mn.nuisance[nidx] + optn.sdev[nidx]*randn()
    if new_nval > optn.bounds[nidx,2] || new_nval < optn.bounds[nidx,1]
        #model should not be modified if prior is violated
        priorviolate = true
    else
        mn.nidx_mem = nidx
        mn.nu_old = mn.nuisance[nidx]
        mn.nuisance[nidx] = new_nval
    end
    statn.move_tries[nidx] += 1
    # @info mn.nu_old, mn.nuisance[nidx], priorviolate
    return 1, priorviolate
end

function undo_move!(mn::ModelNuisance)
    # @info mn.nidx_mem, mn.nu_old, mn.nuisance[mn.nidx_mem]
    mn.nuisance[mn.nidx_mem] = mn.nu_old
    mn.nidx_mem = 0 #delete memory so we can't undo twice
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
            @info typeof(opt), msg
        end
        fill!(stat.move_tries, 0)
        fill!(stat.accepted_moves, 0)
    end
end

function get_acceptance_stats!(isample::Int, opt::OptionsNuisance, stat::Stats)
    if mod(isample-1, opt.stat_window) == 0
        stat.accept_rate[:] = 100. *stat.accepted_moves./stat.move_tries
        if opt.dispstatstoscreen
            @info typeof(opt), stat.accept_rate
        end
        fill!(stat.move_tries, 0)
        fill!(stat.accepted_moves, 0)
    end
end

# history methods
function open_history(opt::Options)
    if isfile(opt.costs_filename)
        @assert (opt.history_mode=="a") "$(opt.costs_filename) exists"
    end
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

#history for nuisance models (which have different internal structure)
function open_history(optn::OptionsNuisance)
    if isfile(optn.costs_filename)
        @assert (optn.history_mode == "a") "$(optn.costs_filename) exists"
    end
    if optn.report_freq > 0
        @info("running nuisance sampler...")
    end
    fp_costs = nothing
    if length(optn.costs_filename) > 0
        fp_costs = open(optn.costs_filename, optn.history_mode)
    end
    fp_vals = nothing
    if length(optn.vals_filename) > 0
        fp_vals = open(optn.vals_filename, optn.history_mode)
    end
    return Writepointers_nuisance(fp_costs, fp_vals)
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

function close_history(wpn::Writepointers_nuisance)
    if wpn.fp_costs != nothing
        close(wpn.fp_costs)
    end
    if wpn.fp_vals != nothing
        close(wpn.fp_vals)
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
    write_history(opt, m.fstar, [m.xtrain; m.ftrain], misfit, stat.accept_rate[1],
                        stat.accept_rate[2], stat.accept_rate[3], stat.accept_rate[4], m.n,
                       isample, wp.fp_costs, wp.fp_fstar, wp.fp_x_ftrain, T, writemodel)
end

function write_history(isample::Int, optn::OptionsNuisance, mn::ModelNuisance, misfit::Float64,
                    statn::Stats, wpn::Writepointers_nuisance, T::Float64, writemodel::Bool)

    write_history(optn, mn.nuisance, misfit, statn.accept_rate, isample, wpn.fp_costs,
                wpn.fp_vals, T, writemodel)
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
        if fp_x_ftrain != nothing
            write(fp_x_ftrain, convert(Array{eltype(Float64)},x_ftrain))
            flush(fp_x_ftrain)
        end
        if writemodel
            if fp_fstar != nothing
                write(fp_fstar, convert(Array{eltype(Float64)},fstar))
                flush(fp_fstar)
            end
        end
    end
end

function write_history(optn::OptionsNuisance, nvals::Array{Float64,1}, misfit::Float64,
                    acceptanceRate::Float64, iter::Int, fp_costs::Union{IOStream, Nothing},
                    fp_vals::Union{IOStream, Nothing}, T::Float64, writemodel::Bool)
    if (mod(iter-1, optn.save_freq) == 0 || iter == 1)
        if fp_costs != nothing
            msg = @sprintf("%d %e %e %e\n", iter, acceptanceRate, misfit, T)
            write(fp_costs, msg)
            flush(fp_costs)
        end
        if fp_vals != nothing
            msg = @sprintf("%d", iter)
            for nval = nvals
                msg *= @sprintf(" %e", nval)
            end
            msg *= "\n"
            write(fp_vals, msg)
            flush(fp_vals)
        end
    end
end

function history(optn::OptionsNuisance; stat=:misfit)
    cols = Dict(:iter => 1, :acceptanceRate => 2,
                :misfit => 3, :T => 4)
    colnum = cols[stat]
    data = readdlm(optn.costs_filename, ' ', Float64)[:,colnum]
    if stat == :iter
        data = Int.(data)
    end
    return data
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
        iters, rem = divrem(filesize(opt.fstar_filename), size(opt.xall,2) * size(opt.fbounds, 1) * sizeof(Float64))
        @assert rem == 0
        fp_models = open(opt.fstar_filename)
        fstar = Array{Array{Float64}}(undef, iters)
        for i = 1:iters
            if typeof(opt) == OptionsStat
                fstar[i] = zeros(Float64, (size(opt.fbounds, 1), size(opt.xall,2)))
            else
                fstar[i] = zeros(Float64, (size(opt.xall,2), size(opt.fbounds, 1)))
            end
            read!(fp_models, fstar[i])
        end
        return fstar
    end
    if stat == :x_ftrain
        if length(opt.x_ftrain_filename) == 0
            @warn("history, requested x_ftrain, but you haven't stored this information.")
            return []
        end
        iters, rem = divrem(filesize(opt.x_ftrain_filename), opt.nmax * sizeof(Float64) * (size(opt.fbounds, 1) +
                                                                                        size(opt.xbounds, 1)))
        @assert rem == 0
        fp_models = open(opt.x_ftrain_filename)
        x_ftrain = Array{Array{Float64,2},1}(undef, iters)
        for i = 1:iters
            x_ftrain[i] = zeros(Float64, size(opt.fbounds, 1) + size(opt.xbounds, 1), opt.nmax)
            read!(fp_models, x_ftrain[i])
        end
        return x_ftrain
    end
    @warn("history, requested stat: $(stat) is not recognized.")
    return []
end
