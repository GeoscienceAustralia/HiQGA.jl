Sys.iswindows() && (ENV["MPLBACKEND"]="qt4agg")
using PyPlot, Test, Random, Revise, Statistics, LinearAlgebra
any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
import GP, TransD_GP
## make options for the model we'll be modifying in McMC
nmin, nmax = 2, 400
λ, δ = [0.1, 0.05], 0.1
fbounds = [-2 2.]
demean = true
sdev_prop = [0.1]
sdev_pos = [0.05;0.05]
pnorm = 2.
K = GP.Mat32()

λx,λy = 0.6,0.6
x = 0:(0.01λx):λx
y = 0:(0.01λy):2λy
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
opt = TransD_GP.Options(nmin = nmin,
                        nmax = nmax,
                        xbounds = xbounds,
                        fbounds = fbounds,
                        xall = xall,
                        λ = λ,
                        δ = δ,
                        demean = demean,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        pnorm = pnorm,
                        quasimultid = false,
                        K = K
                        )
@time m = TransD_GP.init(opt)
## run tests for the different McMC moves
@testset "GP and MCMC move do and undo state tests" begin
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n], opt.xall, opt.λ², opt.δ, nogetvars=true, demean=demean, p=pnorm)
    @testset "init test" begin @test norm(mean(ftest' - 0.5log10.(m.fstar))) < 1e-12 end
    @testset "birth tests" begin
    for i = 1:100
        TransD_GP.birth!(m, opt)
    end
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n], opt.xall, opt.λ², opt.δ, nogetvars=true, demean=demean, p=pnorm)
    @test norm(mean(ftest' - 0.5log10.(m.fstar))) < 1e-12
    end
    @testset "death tests" begin
    for i = 1:100
        TransD_GP.death!(m, opt)
    end
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n], opt.xall, opt.λ², opt.δ, nogetvars=true, demean=demean, p=pnorm)
    @test norm(mean(ftest' - 0.5log10.(m.fstar))) < 1e-12
    end
    # birth and death hold correct states if tests above passed
    @testset "undo birth" begin
    mold = deepcopy(m)
    TransD_GP.birth!(m, opt)
    TransD_GP.undo_birth!(m, opt)
    TransD_GP.sync_model!(m, opt)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n], opt.xall, opt.λ², opt.δ, nogetvars=true, demean=demean, p=pnorm)
    @test norm(mean(ftest' - 0.5log10.(m.fstar))) < 1e-12
    @test norm(mean(mold.fstar - m.fstar)) < 1e-12
    end
    # undo_birth holds state if gotten till here
    @testset "undo multiple birth" begin
    TransD_GP.birth!(m, opt)
    TransD_GP.birth!(m, opt)
    mold = deepcopy(m)
    TransD_GP.birth!(m, opt)
    TransD_GP.undo_birth!(m, opt)
    TransD_GP.sync_model!(m, opt)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n], opt.xall, opt.λ², opt.δ, nogetvars=true, demean=demean, p=pnorm)
    @test norm(mean(ftest' - 0.5log10.(m.fstar))) < 1e-12
    @test norm(mean(mold.fstar - m.fstar)) < 1e-12
    end
    # undo_birth holds state as well for multiple births and deaths till here
    @testset "undo death" begin
    mold = deepcopy(m)
    TransD_GP.death!(m, opt)
    TransD_GP.undo_death!(m, opt)
    TransD_GP.sync_model!(m, opt)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n], opt.xall, opt.λ², opt.δ, nogetvars=true, demean=demean, p=pnorm)
    @test norm(mean(ftest' - 0.5log10.(m.fstar))) < 1e-12
    @test norm(mean(mold.fstar - m.fstar)) < 1e-12
    end
    # undo death holds state if here
    @testset "undo multiple death" begin
    TransD_GP.birth!(m, opt)
    TransD_GP.birth!(m, opt)
    TransD_GP.birth!(m, opt)
    mold = deepcopy(m)
    TransD_GP.death!(m, opt)
    TransD_GP.undo_death!(m, opt)
    TransD_GP.sync_model!(m, opt)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n], opt.xall, opt.λ², opt.δ, nogetvars=true, demean=demean, p=pnorm)
    @test norm(mean(ftest' - 0.5log10.(m.fstar))) < 1e-12
    @test norm(mean(mold.fstar - m.fstar)) < 1e-12
    end
    # undo_death holds state as well for multiple births and deaths till here
    @testset "property change" begin
    TransD_GP.property_change!(m, opt)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n], opt.xall, opt.λ², opt.δ, nogetvars=true, demean=demean, p=pnorm)
    @test norm(mean(ftest' - 0.5log10.(m.fstar))) < 1e-12
    end
    # property change works if here
    @testset "undo property change" begin
    mold = deepcopy(m)
    TransD_GP.property_change!(m, opt)
    TransD_GP.undo_property_change!(m, opt)
    TransD_GP.sync_model!(m, opt)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n], opt.xall, opt.λ², opt.δ, nogetvars=true, demean=demean, p=pnorm)
    @test norm(mean(ftest' - 0.5log10.(m.fstar))) < 1e-12
    @test norm(mean(mold.fstar - m.fstar)) < 1e-12
    end
    # undo property change works if here
    @testset "position change" begin
    TransD_GP.position_change!(m, opt)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n], opt.xall, opt.λ², opt.δ, nogetvars=true, demean=demean, p=pnorm)
    @test norm(mean(ftest' - 0.5log10.(m.fstar))) < 1e-12
    end
    # position change works if here
    @testset "undo position change" begin
    mold = deepcopy(m)
    TransD_GP.position_change!(m, opt)
    TransD_GP.undo_position_change!(m, opt)
    TransD_GP.sync_model!(m, opt)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n], opt.xall, opt.λ², opt.δ, nogetvars=true, demean=demean, p=pnorm)
    @test norm(mean(ftest' - 0.5log10.(m.fstar))) < 1e-12
    @test norm(mean(mold.fstar - m.fstar)) < 1e-12
    end
    # undo position change works if here
end
## timing for birth
@time for i = 1:300
          TransD_GP.birth!(m, opt)
end
##
m.fstar = 0.5log10.(m.fstar)
fig,ax = plt.subplots(1, 1, sharex=true, sharey=true)
vmin, vmax = minimum(m.fstar[1,:]), maximum(m.fstar[1,:])
im1 = ax.imshow(reshape(m.fstar[1,:],length(y), length(x)), extent=[x[1],x[end],y[end],y[1]],
    vmin=vmin, vmax=vmax)
colorbar(im1, ax=ax)
ax.scatter(m.xtrain[1,1:m.n], m.xtrain[2,1:m.n],marker="+",c="r")
ax.scatter(m.xtrain[1,1:m.n], m.xtrain[2,1:m.n],c=m.ftrain[1,1:m.n], alpha=0.8)
## Multichannel Trans-D
Random.seed!(11)
fbounds = [-2 2.; 3 4]
sdev_prop = [0.1, 0.05]
opt = TransD_GP.Options(nmin = nmin,
                        nmax = nmax,
                        xbounds = xbounds,
                        fbounds = fbounds,
                        xall = xall,
                        λ = λ,
                        δ = δ,
                        demean = demean,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        pnorm = pnorm,
                        quasimultid = false,
                        K = K
                        )
@time m = TransD_GP.init(opt)
## run tests for the different McMC moves, multichannel training/test
@testset "MultichannelGP and MCMC move do and undo state tests" begin
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n], opt.xall, opt.λ², opt.δ, nogetvars=true, demean=demean, p=pnorm)
    @testset "init test" begin @test norm(mean(ftest' - 0.5log10.(m.fstar))) < 1e-12 end
    @testset "birth tests" begin
    for i = 1:100
        TransD_GP.birth!(m, opt)
    end
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n], opt.xall, opt.λ², opt.δ, nogetvars=true, demean=demean, p=pnorm)
    @test norm(mean(ftest' - 0.5log10.(m.fstar))) < 1e-12
    end
    @testset "death tests" begin
    for i = 1:100
        TransD_GP.death!(m, opt)
    end
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n], opt.xall, opt.λ², opt.δ, nogetvars=true, demean=demean, p=pnorm)
    @test norm(mean(ftest' - 0.5log10.(m.fstar))) < 1e-12
    end
    # birth and death hold correct states if tests above passed
    @testset "undo birth" begin
    mold = deepcopy(m)
    TransD_GP.birth!(m, opt)
    TransD_GP.undo_birth!(m, opt)
    TransD_GP.sync_model!(m, opt)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n], opt.xall, opt.λ², opt.δ, nogetvars=true, demean=demean, p=pnorm)
    @test norm(mean(ftest' - 0.5log10.(m.fstar))) < 1e-12
    @test norm(mean(mold.fstar - m.fstar)) < 1e-12
    end
    # undo_birth holds state if gotten till here
    @testset "undo multiple birth" begin
    TransD_GP.birth!(m, opt)
    TransD_GP.birth!(m, opt)
    mold = deepcopy(m)
    TransD_GP.birth!(m, opt)
    TransD_GP.undo_birth!(m, opt)
    TransD_GP.sync_model!(m, opt)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n], opt.xall, opt.λ², opt.δ, nogetvars=true, demean=demean, p=pnorm)
    @test norm(mean(ftest' - 0.5log10.(m.fstar))) < 1e-12
    @test norm(mean(mold.fstar - m.fstar)) < 1e-12
    end
    # undo_birth holds state as well for multiple births and deaths till here
    @testset "undo death" begin
    mold = deepcopy(m)
    TransD_GP.death!(m, opt)
    TransD_GP.undo_death!(m, opt)
    TransD_GP.sync_model!(m, opt)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n], opt.xall, opt.λ², opt.δ, nogetvars=true, demean=demean, p=pnorm)
    @test norm(mean(ftest' - 0.5log10.(m.fstar))) < 1e-12
    @test norm(mean(mold.fstar - m.fstar)) < 1e-12
    end
    # undo death holds state if here
    @testset "undo multiple death" begin
    TransD_GP.birth!(m, opt)
    TransD_GP.birth!(m, opt)
    TransD_GP.birth!(m, opt)
    mold = deepcopy(m)
    TransD_GP.death!(m, opt)
    TransD_GP.undo_death!(m, opt)
    TransD_GP.sync_model!(m, opt)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n], opt.xall, opt.λ², opt.δ, nogetvars=true, demean=demean, p=pnorm)
    @test norm(mean(ftest' - 0.5log10.(m.fstar))) < 1e-12
    @test norm(mean(mold.fstar - m.fstar)) < 1e-12
    end
    # undo_death holds state as well for multiple births and deaths till here
    @testset "property change" begin
    TransD_GP.property_change!(m, opt)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n], opt.xall, opt.λ², opt.δ, nogetvars=true, demean=demean, p=pnorm)
    @test norm(mean(ftest' - 0.5log10.(m.fstar))) < 1e-12
    end
    # property change works if here
    @testset "undo property change" begin
    mold = deepcopy(m)
    TransD_GP.property_change!(m, opt)
    TransD_GP.undo_property_change!(m, opt)
    TransD_GP.sync_model!(m, opt)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n], opt.xall, opt.λ², opt.δ, nogetvars=true, demean=demean, p=pnorm)
    @test norm(mean(ftest' - 0.5log10.(m.fstar))) < 1e-12
    @test norm(mean(mold.fstar - m.fstar)) < 1e-12
    end
    # undo property change works if here
    @testset "position change" begin
    TransD_GP.position_change!(m, opt)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n], opt.xall, opt.λ², opt.δ, nogetvars=true, demean=demean, p=pnorm)
    @test norm(mean(ftest' - 0.5log10.(m.fstar))) < 1e-12
    end
    # position change works if here
    @testset "undo position change" begin
    mold = deepcopy(m)
    TransD_GP.position_change!(m, opt)
    TransD_GP.undo_position_change!(m, opt)
    TransD_GP.sync_model!(m, opt)
    ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n], opt.xall, opt.λ², opt.δ, nogetvars=true, demean=demean, p=pnorm)
    @test norm(mean(ftest' - 0.5log10.(m.fstar))) < 1e-12
    @test norm(mean(mold.fstar - m.fstar)) < 1e-12
    end
    # undo position change works if here
end
## timing for birth
    @time for i = 1:300
              TransD_GP.birth!(m, opt)
    end
## plot
    m.fstar = 0.5log10.(m.fstar)
    fig,ax = plt.subplots(1,2, sharex=true, sharey=true)
    vmin, vmax = minimum(m.fstar[1,:]), maximum(m.fstar[1,:])
    im1 = ax[1].imshow(reshape(m.fstar[1,:],length(y), length(x)), extent=[x[1],x[end],y[end],y[1]],
        vmin=vmin, vmax=vmax)
    colorbar(im1, ax=ax[1])
    ax[1].scatter(m.xtrain[1,1:m.n], m.xtrain[2,1:m.n],marker="+",c="r")
    ax[1].scatter(m.xtrain[1,1:m.n], m.xtrain[2,1:m.n],c=m.ftrain[1,1:m.n], alpha=0.8)
    vmin, vmax = minimum(m.fstar[2,:]), maximum(m.fstar[2,:])
    im2 = ax[2].imshow(reshape(m.fstar[2,:],length(y), length(x)), extent=[x[1],x[end],y[end],y[1]],
        vmin=vmin, vmax=vmax)
    colorbar(im2, ax=ax[2])
    ax[2].scatter(m.xtrain[1,1:m.n], m.xtrain[2,1:m.n],marker="+",c="r")
    ax[2].scatter(m.xtrain[1,1:m.n], m.xtrain[2,1:m.n],c=m.ftrain[2,1:m.n], alpha=0.8,
        vmin=vmin, vmax=vmax)
## but we also need to check if the μ model is restored
