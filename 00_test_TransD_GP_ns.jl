Sys.iswindows() && (ENV["MPLBACKEND"]="qt4agg")
using PyPlot, Test, Random, Revise, Statistics, LinearAlgebra
any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
import GP, TransD_GP
## make options for the multichannel lengthscale GP
nminlog10λ, nmaxlog10λ = 2, 100
pnorm = 2.
Klog10λ = GP.Mat32()
λx,λy = 0.6,0.6
x = 0:(0.01λx):λx
y = 0:(0.01λy):2λy
λlog10λ = [0.1maximum(y), 0.1maximum(x)]
demean = true
sdev_poslog10λ = [0.01maximum(y), 0.01maximum(x)]
log10bounds = [-2 -0.69; -2 -0.69]
δlog10λ = 0.1
sdev_proplog10λ = [0.1, 0.1]
xall = zeros(2,length(x)*length(y))
for i in 1:size(xall,2)
    yid, xid = Tuple(CartesianIndices((length(y),length(x)))[i])
    xall[:,i] = [x[xid]; y[yid]]
end
xbounds = zeros(Float64,size(xall, 1), 2)
for dim in 1:size(xall, 1)
    xbounds[dim,:] = [minimum(xall[dim,:]), maximum(xall[dim,:])]
end
## Initialize a model using these options
Random.seed!(12)
optlog10λ = TransD_GP.Options(nmin = nminlog10λ,
                        nmax = nmaxlog10λ,
                        xbounds = xbounds,
                        fbounds = log10bounds,
                        xall = xall,
                        λ = λlog10λ,
                        δ = δlog10λ,
                        demean = demean,
                        sdev_prop = sdev_proplog10λ,
                        sdev_pos = sdev_poslog10λ,
                        pnorm = pnorm,
                        quasimultid = false,
                        K = Klog10λ
                        )
@time  log10λ = TransD_GP.init(optlog10λ)
for i = 1:49
    TransD_GP.birth!(log10λ, optlog10λ)
end
## make options for the nonstationary GP
nmin, nmax = 2, 400
fbounds = [-2. 2]
δ = 0.1
sdev_prop = [0.1]
sdev_pos = [0.05, 0.05]
K = GP.Mat32()
## Initialize model for the nonstationary GP
opt = TransD_GP.Options(nmin = nmin,
                        nmax = nmax,
                        xbounds = optlog10λ.xbounds,
                        fbounds = fbounds,
                        xall = optlog10λ.xall,
                        λ = λlog10λ,
                        δ = δ,
                        demean = demean,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        pnorm = pnorm,
                        quasimultid = false,
                        K = K
                        )
## Nonstationary GP model changes
λ² = log10λ.fstar
@time m = TransD_GP.init(opt, log10λ)
@testset "Nonstationary GP model changes" begin
    @testset "init test" begin
    idxs = TransD_GP.gettrainidx(opt.kdtree, m.xtrain, m.n)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n],
        opt.xall, λ², λ²[:,idxs], δ, p=2, demean=demean, nogetvars=true)
    @test norm(mean(ftest - m.fstar)) < 1e-12
    end
    @testset "birth tests" begin
    for i = 1:100
        TransD_GP.birth!(m, opt, log10λ)
    end
    idxs = TransD_GP.gettrainidx(opt.kdtree, m.xtrain, m.n)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n],
        opt.xall, λ², λ²[:,idxs], δ, p=2, demean=demean, nogetvars=true)
    @test norm(mean(ftest - m.fstar)) < 1e-12
    end
    @testset "death tests" begin
    for i = 1:100
        TransD_GP.death!(m, opt)
    end
    idxs = TransD_GP.gettrainidx(opt.kdtree, m.xtrain, m.n)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n],
        opt.xall, λ², λ²[:,idxs], δ, p=2, demean=demean, nogetvars=true)
    @test norm(mean(ftest - m.fstar)) < 1e-12
    end
    @testset "undo birth" begin
    mold = deepcopy(m)
    TransD_GP.birth!(m, opt, log10λ)
    TransD_GP.undo_birth!(m, opt)
    TransD_GP.sync_model!(m, opt)
    idxs = TransD_GP.gettrainidx(opt.kdtree, m.xtrain, m.n)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n],
        opt.xall, λ², λ²[:,idxs], δ, p=2, demean=demean, nogetvars=true)
    @test norm(mean(ftest - m.fstar)) < 1e-12
    @test norm(mean(mold.fstar - m.fstar)) < 1e-12
    end
    @testset "undo multiple birth" begin
    TransD_GP.birth!(m, opt, log10λ )
    TransD_GP.birth!(m, opt, log10λ)
    mold = deepcopy(m)
    TransD_GP.birth!(m, opt, log10λ)
    TransD_GP.undo_birth!(m, opt)
    TransD_GP.sync_model!(m, opt)
    idxs = TransD_GP.gettrainidx(opt.kdtree, m.xtrain, m.n)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n],
        opt.xall, λ², λ²[:,idxs], δ, p=2, demean=demean, nogetvars=true)
    @test norm(mean(ftest - m.fstar)) < 1e-12
    @test norm(mean(mold.fstar - m.fstar)) < 1e-12
    end
    @testset "undo death" begin
    mold = deepcopy(m)
    TransD_GP.death!(m, opt)
    TransD_GP.undo_death!(m, opt)
    TransD_GP.sync_model!(m, opt)
    idxs = TransD_GP.gettrainidx(opt.kdtree, m.xtrain, m.n)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n],
        opt.xall, λ², λ²[:,idxs], δ, p=2, demean=demean, nogetvars=true)
    @test norm(mean(ftest - m.fstar)) < 1e-12
    @test norm(mean(mold.fstar - m.fstar)) < 1e-12
    end
    @testset "undo multiple death" begin
    TransD_GP.birth!(m, opt, log10λ)
    TransD_GP.birth!(m, opt, log10λ)
    TransD_GP.birth!(m, opt, log10λ)
    mold = deepcopy(m)
    TransD_GP.death!(m, opt)
    TransD_GP.undo_death!(m, opt)
    TransD_GP.sync_model!(m, opt)
    idxs = TransD_GP.gettrainidx(opt.kdtree, m.xtrain, m.n)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n],
        opt.xall, λ², λ²[:,idxs], δ, p=2, demean=demean, nogetvars=true)
    @test norm(mean(ftest - m.fstar)) < 1e-12
    @test norm(mean(mold.fstar - m.fstar)) < 1e-12
    end
    @testset "property change" begin
    TransD_GP.property_change!(m, opt)
    idxs = TransD_GP.gettrainidx(opt.kdtree, m.xtrain, m.n)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n],
        opt.xall, λ², λ²[:,idxs], δ, p=2, demean=demean, nogetvars=true)
    @test norm(mean(ftest - m.fstar)) < 1e-12
    end
    @testset "undo property change" begin
    mold = deepcopy(m)
    TransD_GP.property_change!(m, opt)
    TransD_GP.undo_property_change!(m, opt)
    TransD_GP.sync_model!(m, opt)
    idxs = TransD_GP.gettrainidx(opt.kdtree, m.xtrain, m.n)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n],
        opt.xall, λ², λ²[:,idxs], δ, p=2, demean=demean, nogetvars=true)
    @test norm(mean(ftest - m.fstar)) < 1e-12
    @test norm(mean(mold.fstar - m.fstar)) < 1e-12
    end
    @testset "position change" begin
    TransD_GP.position_change!(m, opt, log10λ)
    idxs = TransD_GP.gettrainidx(opt.kdtree, m.xtrain, m.n)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n],
        opt.xall, λ², λ²[:,idxs], δ, p=2, demean=demean, nogetvars=true)
    @test norm(mean(ftest - m.fstar)) < 1e-12
    end
    @testset "undo position change" begin
    mold = deepcopy(m)
    TransD_GP.position_change!(m, opt, log10λ)
    TransD_GP.undo_position_change!(m, opt, log10λ)
    TransD_GP.sync_model!(m, opt)
    idxs = TransD_GP.gettrainidx(opt.kdtree, m.xtrain, m.n)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n],
        opt.xall, λ², λ²[:,idxs], δ, p=2, demean=demean, nogetvars=true)
    @test norm(mean(ftest - m.fstar)) < 1e-12
    @test norm(mean(mold.fstar - m.fstar)) < 1e-12
    end
end
## timing for birth
@time for i = 1:300
    TransD_GP.birth!(m, opt, log10λ)
end
## plot
fig,ax = plt.subplots(1,3, sharex=true, sharey=true)
log10λ.fstar = sqrt.(log10λ.fstar)
vmin, vmax = minimum(log10λ.fstar[1,:]), maximum(log10λ.fstar[1,:])
im1 = ax[1].imshow(reshape(log10λ.fstar[1,:],length(y), length(x)), extent=[x[1],x[end],y[end],y[1]],
    vmin=vmin, vmax=vmax)
colorbar(im1, ax=ax[1])
ax[1].scatter(log10λ.xtrain[1,1:log10λ.n], log10λ.xtrain[2,1:log10λ.n],marker="+",c="r")
ax[1].scatter(log10λ.xtrain[1,1:log10λ.n], log10λ.xtrain[2,1:log10λ.n],c=log10λ.ftrain[1,1:log10λ.n], alpha=0.8)
vmin, vmax = minimum(log10λ.fstar[2,:]), maximum(log10λ.fstar[2,:])
im2 = ax[2].imshow(reshape(log10λ.fstar[2,:],length(y), length(x)), extent=[x[1],x[end],y[end],y[1]],
    vmin=vmin, vmax=vmax)
colorbar(im2, ax=ax[2])
ax[2].scatter(log10λ.xtrain[1,1:log10λ.n], log10λ.xtrain[2,1:log10λ.n],marker="+",c="r")
ax[2].scatter(log10λ.xtrain[1,1:log10λ.n], log10λ.xtrain[2,1:log10λ.n],c=log10λ.ftrain[2,1:log10λ.n], alpha=0.8,
    vmin=vmin, vmax=vmax)
vmin, vmax = minimum(m.fstar), maximum(m.fstar)
im3 = ax[3].imshow(reshape(m.fstar,length(y), length(x)), extent=[x[1],x[end],y[end],y[1]],
    vmin=vmin, vmax=vmax)
colorbar(im3, ax=ax[3])
ax[3].scatter(m.xtrain[1,1:m.n], m.xtrain[2,1:m.n],marker="+",c="r")
ax[3].scatter(m.xtrain[1,1:m.n], m.xtrain[2,1:m.n],c=m.ftrain[1:m.n], alpha=0.8,
        vmin=vmin, vmax=vmax)
