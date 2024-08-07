## make options for the stationary multichannel lengthscale GP
nminlog10λ, nmaxlog10λ = 2, 200
pnorm = 2.
Klog10λ = transD_GP.GP.Mat32()
λx,λy = 1, 1
x = 0:(0.005λx):λx
y = 0:(0.005λy):0.5λy
λlog10λ = [0.05maximum(x), 0.05maximum(y)]
sdev_poslog10λ = [0.01maximum(y), 0.01maximum(x)]
log10bounds = [-1 -0.69; -1 -0.69]
δlog10λ = 0.1
sdev_proplog10λ = [0.1, 0.1]
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
optlog10λ = transD_GP.OptionsStat(nmin = nminlog10λ,
                        nmax = nmaxlog10λ,
                        xbounds = xbounds,
                        fbounds = log10bounds,
                        xall = xall,
                        λ = λlog10λ,
                        δ = δlog10λ,
                        sdev_prop = sdev_proplog10λ,
                        sdev_pos = sdev_poslog10λ,
                        pnorm = pnorm,
                        quasimultid = false,
                        K = Klog10λ,
                        needλ²fromlog = true,
                        updatenonstat = true,
                        sampledc = false, # cannot sampledc for length scales
                        timesλ = 3.6
                        )
@time  log10λ = transD_GP.init(optlog10λ, 1)
## make options for the nonstationary GP -- changes to these should not 
# affect the stationary GP
nmin, nmax = 2, 200
fbounds = [-2. 2]
δ = 0.1
sdev_prop = [0.1]
sdev_pos = [0.05, 0.05]
K = transD_GP.GP.Mat32()
demean_ns = false
sampledc_ns = true
## Initialize model for the nonstationary GP
Random.seed!(13)
opt = transD_GP.OptionsNonstat(optlog10λ,
                        nmin = nmin,
                        nmax = nmax,
                        fbounds = fbounds,
                        δ = δ,
                        demean = demean_ns,
                        sampledc = sampledc_ns,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        pnorm = pnorm,
                        K = K
                        )
@time m = transD_GP.init(opt, log10λ, 1)
##
fracthresh = 0.05
lsc = transD_GP.testupdate(optlog10λ, log10λ)
ftest = transD_GP.testupdate(opt, log10λ, m)
@testset "Nonstationary GP model changes" begin
    @testset "init test" begin
        @test norm(mean(lsc' - 0.5log10.(log10λ.fstar))) < 1e-12
        for i = 1:150
            transD_GP.birth!(log10λ, optlog10λ, m, opt)
        end
        lsc[:] = transD_GP.testupdate(optlog10λ, log10λ)
        ftest = transD_GP.testupdate(opt, log10λ, m)
        @test norm(mean(lsc' - 0.5log10.(log10λ.fstar))) < 1e-12
        @test mean(abs.((ftest - m.fstar)./ftest)) < fracthresh
    end
    @testset "birth tests" begin
        for i = 1:100
            transD_GP.birth!(m, opt, log10λ)
        end
        ftest = transD_GP.testupdate(opt, log10λ, m)
        @test norm(mean(lsc' - 0.5log10.(log10λ.fstar))) < 1e-12
        @test mean(abs.((ftest - m.fstar)./ftest)) < fracthresh
    end
    @testset "death tests" begin
        for i = 1:100
            transD_GP.death!(m, opt)
        end
        ftest = transD_GP.testupdate(opt, log10λ, m)
        @test norm(mean(lsc' - 0.5log10.(log10λ.fstar))) < 1e-12
        @test mean(abs.((ftest - m.fstar)./ftest)) < fracthresh
    end
    @testset "undo birth" begin
        mold = deepcopy(m)
        transD_GP.birth!(m, opt, log10λ)
        transD_GP.undo_birth!(m)
        transD_GP.sync_model!(m, opt)
        ftest = transD_GP.testupdate(opt, log10λ, m)
        @test norm(mean(lsc' - 0.5log10.(log10λ.fstar))) < 1e-12
        @test mean(abs.((ftest - m.fstar)./ftest)) < fracthresh
        @test norm(mean(mold.fstar - m.fstar)) < 1e-12
    end
    @testset "undo multiple birth" begin
        transD_GP.birth!(m, opt, log10λ )
        transD_GP.birth!(m, opt, log10λ)
        mold = deepcopy(m)
        transD_GP.birth!(m, opt, log10λ)
        transD_GP.undo_birth!(m)
        transD_GP.sync_model!(m, opt)
        ftest = transD_GP.testupdate(opt, log10λ, m)
        @test norm(mean(lsc' - 0.5log10.(log10λ.fstar))) < 1e-12
        @test mean(abs.((ftest - m.fstar)./ftest)) < fracthresh
        @test norm(mean(mold.fstar - m.fstar)) < 1e-12
    end
    @testset "undo death" begin
        mold = deepcopy(m)
        transD_GP.death!(m, opt)
        transD_GP.undo_death!(m, opt)
        transD_GP.sync_model!(m, opt)
        ftest = transD_GP.testupdate(opt, log10λ, m)
        @test norm(mean(lsc' - 0.5log10.(log10λ.fstar))) < 1e-12
        @test mean(abs.((ftest - m.fstar)./ftest)) < fracthresh
        @test norm(mean(mold.fstar - m.fstar)) < 1e-12
    end
    @testset "undo multiple death" begin
        transD_GP.birth!(m, opt, log10λ)
        transD_GP.birth!(m, opt, log10λ)
        transD_GP.birth!(m, opt, log10λ)
        mold = deepcopy(m)
        transD_GP.death!(m, opt)
        transD_GP.undo_death!(m, opt)
        transD_GP.sync_model!(m, opt)
        ftest = transD_GP.testupdate(opt, log10λ, m)
        @test norm(mean(lsc' - 0.5log10.(log10λ.fstar))) < 1e-12
        @test mean(abs.((ftest - m.fstar)./ftest)) < fracthresh
        @test norm(mean(mold.fstar - m.fstar)) < 1e-12
    end
    @testset "property change" begin
        transD_GP.property_change!(m, opt)
        ftest = transD_GP.testupdate(opt, log10λ, m)
        @test norm(mean(lsc' - 0.5log10.(log10λ.fstar))) < 1e-12
        @test mean(abs.((ftest - m.fstar)./ftest)) < fracthresh
    end
    @testset "undo property change" begin
        mold = deepcopy(m)
        transD_GP.property_change!(m, opt)
        transD_GP.undo_property_change!(m)
        transD_GP.sync_model!(m, opt)
        ftest = transD_GP.testupdate(opt, log10λ, m)
        @test norm(mean(lsc' - 0.5log10.(log10λ.fstar))) < 1e-12
        @test mean(abs.((ftest - m.fstar)./ftest)) < fracthresh
        @test norm(mean(mold.fstar - m.fstar)) < 1e-12
    end
    @testset "position change" begin
        transD_GP.position_change!(m, opt, log10λ)
        ftest = transD_GP.testupdate(opt, log10λ, m)
        @test norm(mean(lsc' - 0.5log10.(log10λ.fstar))) < 1e-12
        @test mean(abs.((ftest - m.fstar)./ftest)) < fracthresh
    end
    @testset "undo position change" begin
        mold = deepcopy(m)
        transD_GP.position_change!(m, opt, log10λ)
        transD_GP.undo_position_change!(m, opt, log10λ)
        transD_GP.sync_model!(m, opt)
        ftest = transD_GP.testupdate(opt, log10λ, m)
        @test norm(mean(lsc' - 0.5log10.(log10λ.fstar))) < 1e-12
        @test mean(abs.((ftest - m.fstar)./ftest)) < fracthresh
        @test norm(mean(mold.fstar - m.fstar)) < 1e-12
    end
    @testset "dc change" begin
        transD_GP.dc_change!(m, opt)
        transD_GP.sync_model!(m, opt)
        ftest = transD_GP.testupdate(opt, log10λ, m)
        @test norm(mean(lsc' - 0.5log10.(log10λ.fstar))) < 1e-12
        @test mean(abs.((ftest - m.fstar)./ftest)) < fracthresh
    end
    @testset "undo dc change" begin
        mold = deepcopy(m)
        transD_GP.dc_change!(m, opt)
        transD_GP.undo_dc_change!(m)
        transD_GP.sync_model!(m, opt)
        ftest = transD_GP.testupdate(opt, log10λ, m)
        @test norm(mean(lsc' - 0.5log10.(log10λ.fstar))) < 1e-12
        @test mean(abs.((ftest - m.fstar)./ftest)) < fracthresh
        @test norm(mean(mold.fstar - m.fstar)) < 1e-12
    end
end
# ## timing for ns birth in ns model
# log10λ = transD_GP.init(optlog10λ)
# m = transD_GP.init(opt, log10λ)
# for i = 1:148
#     transD_GP.birth!(log10λ, optlog10λ, m, opt)
# end
# @time for i = 1:148
#     transD_GP.birth!(m, opt, log10λ)
# end
# ## time for ntimes birth death
# NTIMES = 20
# ntimes = 150
# T = zeros(NTIMES)
# for I = 1:NTIMES
#     T[I] = time()
#     for i = 1:ntimes
#         transD_GP.birth!(m, opt, log10λ)
#         transD_GP.death!(m, opt)
#     end
#     T[I] = (time() - T[I])/2ntimes
# end
# @info "time for $ntimes birth/death is $(mean(T)) +- $(std(T)/sqrt(NTIMES))"
# ## plot
# l = log10.(log10λ.fstar.^0.5)
# fig,ax = plt.subplots(1,3, sharex=true, sharey=true)
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
# im3 = ax[3].imshow(reshape(m.fstar,length(x), length(y)), extent=[y[1],y[end],x[end],x[1]],
#     vmin=vmin, vmax=vmax)
# colorbar(im3, ax=ax[3])
# ax[3].scatter(m.xtrain[2,1:m.n], m.xtrain[1,1:m.n], marker="+",c="r")
# ax[3].scatter(m.xtrain[2,1:m.n], m.xtrain[1,1:m.n], c=m.ftrain[1:m.n], alpha=0.8,
#         vmin=vmin, vmax=vmax)
