## make options for a purely stationary properties GP
nmin, nmax = 2, 200
pnorm = 2.
λx, λy = 1, 1
x = 0:(0.005λx):λx
y = 0:(0.005λy):0.5λy
nmin, nmax = 2, 200
fbounds = [-2. 2]
δ = 0.1
sdev_prop = [0.1]
sdev_pos = [0.05, 0.05]
K = transD_GP.GP.Mat32()
λ = [0.05maximum(x), 0.05maximum(y)]
demean = false
sampledc = true
##
xall = zeros(2,length(x)*length(y))
for i in 1:size(xall,2)
    xid, yid = Tuple(CartesianIndices((length(x),length(y)))[i])
    xall[:,i] = [x[xid]; y[yid]]
end
xbounds = zeros(Float64,size(xall, 1), 2)
for dim in 1:size(xall, 1)
    xbounds[dim,:] = [minimum(xall[dim,:]), maximum(xall[dim,:])]
end
## Initialize a model using these options
Random.seed!(12)
opt = transD_GP.OptionsStat(nmin = nmin,
                        nmax = nmax,
                        xbounds = xbounds,
                        fbounds = fbounds,
                        xall = xall,
                        λ = λ,
                        δ = δ,
                        demean = demean,
                        sampledc = sampledc,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        pnorm = pnorm,
                        needλ²fromlog = false,
                        updatenonstat = false,
                        quasimultid = false,
                        K = K,
                        timesλ = 3,
                        )
@time m = transD_GP.init(opt)
## run tests for the different McMC moves
@testset "Stationary GP and MCMC move do and undo state tests" begin
    @testset "init test" begin
        ftest = transD_GP.testupdate(opt, m)
        @test norm(mean(ftest' - m.fstar)) < 1e-12
    end
    @testset "birth tests" begin
        for i = 1:50
            transD_GP.birth!(m, opt)
        end
        ftest = transD_GP.testupdate(opt, m)
        @test norm(mean(ftest' - m.fstar)) < 1e-12
    end
    @testset "death tests" begin
        for i = 1:20
            transD_GP.death!(m, opt)
        end
        ftest = transD_GP.testupdate(opt, m)
        @test norm(mean(ftest' - m.fstar)) < 1e-12
    end
    # birth and death hold correct states if tests above passed
    @testset "undo birth" begin
        transD_GP.birth!(m, opt)
        transD_GP.undo_birth!(m)
        transD_GP.sync_model!(m, opt)
        ftest = transD_GP.testupdate(opt, m)
        @test norm(mean(ftest' - m.fstar)) < 1e-12
    end
    # undo_birth holds state if gotten till here
    @testset "undo multiple birth" begin
        transD_GP.birth!(m, opt)
        transD_GP.birth!(m, opt)
        transD_GP.birth!(m, opt)
        transD_GP.undo_birth!(m)
        transD_GP.sync_model!(m, opt)
        ftest = transD_GP.testupdate(opt, m)
        @test norm(mean(ftest' - m.fstar)) < 1e-12
    end
    # undo_birth holds state as well for multiple births and deaths till here
    @testset "undo death" begin
        transD_GP.death!(m, opt)
        transD_GP.undo_death!(m, opt)
        transD_GP.sync_model!(m, opt)
        ftest = transD_GP.testupdate(opt, m)
        @test norm(mean(ftest' - m.fstar)) < 1e-12
    end
    # undo death holds state if here
    @testset "undo multiple death" begin
        transD_GP.birth!(m, opt)
        transD_GP.birth!(m, opt)
        transD_GP.birth!(m, opt)
        transD_GP.death!(m, opt)
        transD_GP.death!(m, opt)
        transD_GP.undo_death!(m, opt)
        transD_GP.sync_model!(m, opt)
        ftest = transD_GP.testupdate(opt, m)
        @test norm(mean(ftest' - m.fstar)) < 1e-12
    end
    # undo_death holds state as well for multiple births and deaths till here
    @testset "property change" begin
        transD_GP.property_change!(m, opt)
        ftest = transD_GP.testupdate(opt, m)
        @test norm(mean(ftest' - m.fstar)) < 1e-12
    end
    # property change works if here
    @testset "undo property change" begin
        transD_GP.property_change!(m, opt)
        transD_GP.undo_property_change!(m)
        transD_GP.sync_model!(m, opt)
        ftest = transD_GP.testupdate(opt, m)
        @test norm(mean(ftest' - m.fstar)) < 1e-12
    end
    # undo property change works if here
    @testset "position change" begin
        transD_GP.position_change!(m, opt)
        ftest = transD_GP.testupdate(opt, m)
        @test norm(mean(ftest' - m.fstar)) < 1e-12
    end
    # position change works if here
    @testset "undo position change" begin
        transD_GP.position_change!(m, opt)
        transD_GP.undo_position_change!(m, opt)
        transD_GP.sync_model!(m, opt)
        ftest = transD_GP.testupdate(opt, m)
        @test norm(mean(ftest' - m.fstar)) < 1e-12
    end
    # undo position change works if here
    @testset "dc change" begin
        transD_GP.dc_change!(m, opt)
        transD_GP.sync_model!(m, opt)
        ftest = transD_GP.testupdate(opt, m)
        @test norm(mean(ftest' - m.fstar)) < 1e-12
    end
    @testset "undo dc change" begin
        transD_GP.dc_change!(m, opt)
        transD_GP.undo_dc_change!(m)
        transD_GP.sync_model!(m, opt)
        ftest = transD_GP.testupdate(opt, m)
        @test norm(mean(ftest' - m.fstar)) < 1e-12
    end
end
# ## timing for birth
# m = transD_GP.init(opt, log10λ)
# log10λ = transD_GP.init(optlog10λ)
# @time for i = 1:148
#           transD_GP.birth!(log10λ, optlog10λ, m, opt)
# end
# ## timing for many births
# NTIMES = 20
# ntimes = 150
# T = zeros(NTIMES)
# for I = 1:NTIMES
#     T[I] = time()
#     for i = 1:ntimes
#         transD_GP.death!(log10λ, optlog10λ, m, opt)
#         transD_GP.birth!(log10λ, optlog10λ, m, opt)
#     end
#     T[I] = (time() - T[I])/2ntimes
# end
# @info "time for $ntimes birth/death is $(mean(T)) +- $(std(T)/sqrt(NTIMES))"
# ## plot
# l = m.fstar
# fig,ax = plt.subplots(1,2, sharex=true, sharey=true)
# vmin, vmax = minimum(l[1,:]), maximum(l[1,:])
# im1 = ax[1].imshow(reshape(l[1,:],length(x), length(y)), extent=[y[1],y[end],x[end],x[1]],
#     vmin=vmin, vmax=vmax)
# colorbar(im1, ax=ax[1])
# ax[1].scatter(log10λ.xtrain[2,1:log10λ.n], log10λ.xtrain[1,1:log10λ.n], marker="+",c="r")
# ax[1].scatter(log10λ.xtrain[2,1:log10λ.n], log10λ.xtrain[1,1:log10λ.n], c=log10λ.ftrain[1,1:log10λ.n], alpha=0.8)
# vmin, vmax = minimum(l[2,:]), maximum(l[2,:])
# im2 = ax[2].imshow(reshape(l[2,:],length(x), length(y)), extent=[y[1],y[end],x[end],x[1]],
#     vmin=vmin, vmax=vmax)
# colorbar(im2, ax=ax[2])
# ax[2].scatter(log10λ.xtrain[2,1:log10λ.n], log10λ.xtrain[1,1:log10λ.n], marker="+",c="r")
# ax[2].scatter(log10λ.xtrain[2,1:log10λ.n], log10λ.xtrain[1,1:log10λ.n], c=log10λ.ftrain[2,1:log10λ.n], alpha=0.8,
#     vmin=vmin, vmax=vmax)
# vmin, vmax = minimum(m.fstar), maximum(m.fstar)