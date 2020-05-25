any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
using PyPlot, TransD_GP, Distributed, MCMC_Driver, GP, Random, Revise, GeophysOperator
##
img =     Img(
          filename         = "4.2.01.png",
          dx               = 10.0,
          fractrain        = 0.02,
          dec              = 2,
          gausskernelwidth = 7)
##
img.σ, ftrain, Xtrain =  get_training_data(img,
                   sdmaxfrac = 0.05,
                   ybreak = 1000,
                   takeevery = 4)

plot_data(ftrain, Xtrain, img)

xall = get_all_prediction_points(img)
## stationary GP options
nmin, nmax= 2, 100
pnorm = 2.
K = GP.Mat32()
fbounds = permutedims([extrema(ftrain)...])
δ = 0.3
λ = [141,141]
demean = true
sdev_prop = 0.05*diff(fbounds, dims=2)[:]
sdev_pos = 0.05*abs.([diff([extrema(xall[1,:])...])[1], diff([extrema(xall[2,:])...])[1]])
updatenonstat = false
needλ²fromlog = false
## Initialize a stationary GP using these options
Random.seed!(12)
opt = TransD_GP.OptionsStat(nmin = nmin,
                        nmax = nmax,
                        xbounds = [img.x[1] img.x[end];img.y[1] img.y[end]],
                        fbounds = fbounds,
                        xall = xall,
                        λ = λ,
                        δ = δ,
                        demean = demean,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        pnorm = pnorm,
                        quasimultid = false,
                        K = K,
                        timesλ = 3,
                        needλ²fromlog = needλ²fromlog,
                        updatenonstat = updatenonstat
                        )
## Initialize model for the dummy nonstationary GP
Random.seed!(13)
optdummy = TransD_GP.OptionsNonstat(opt,
                        nmin = nmin,
                        nmax = nmax,
                        fbounds = fbounds,
                        δ = δ,
                        demean = demean,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        pnorm = pnorm,
                        K = K
                        )
calc_simple_RMS(img)

## set up McMC
nsamples, nchains, nchainsatone = 50001, 8, 1
Tmax = 2.50
addprocs(nchains)
@info "workers are $(workers())"
@everywhere any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
@everywhere using Distributed
@everywhere import MCMC_Driver
## run McMCopt.his
@time MCMC_Driver.main(opt, optdummy, img, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
rmprocs(workers())
## plot
GeophysOperator.getchi2forall(opt)
