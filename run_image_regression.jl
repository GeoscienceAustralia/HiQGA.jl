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
                   takeevery = 8)

plot_data(ftrain, Xtrain, img)

xall = get_all_prediction_points(img)
## lengthscale GP options
nminlog10λ, nmaxlog10λ = 2, 100
pnorm = 2.
Klog10λ = GP.Mat32()
log10bounds = [log10(10) log10(2000); log10(10) log10(2000)]
δlog10λ = 0.3
λlog10λ = 0.1*abs.([diff([extrema(xall[1,:])...])[1], diff([extrema(xall[2,:])...])[1]])
demean = false
sdev_proplog10λ = 0.05*diff(log10bounds, dims=2)[:]
sdev_poslog10λ = 0.05*abs.([diff([extrema(xall[1,:])...])[1], diff([extrema(xall[2,:])...])[1]])
## Initialize a lengthscale model using these options
Random.seed!(12)
optlog10λ = TransD_GP.OptionsStat(nmin = nminlog10λ,
                        nmax = nmaxlog10λ,
                        xbounds = [img.x[1] img.x[end];img.y[1] img.y[end]],
                        fbounds = log10bounds,
                        xall = xall,
                        λ = λlog10λ,
                        δ = δlog10λ,
                        demean = demean,
                        sdev_prop = sdev_proplog10λ,
                        sdev_pos = sdev_poslog10λ,
                        pnorm = pnorm,
                        quasimultid = false,
                        K = Klog10λ,
                        timesλ = 3
                        )
## make options for the nonstationary actual properties GP
nmin, nmax = 2, 100
fbounds = permutedims([extrema(ftrain)...])
sdev_prop = 0.05*diff(fbounds, dims=2)[:]
sdev_pos = 0.05*abs.([diff([extrema(xall[1,:])...])[1], diff([extrema(xall[2,:])...])[1]])
demean_ns = true
δ = 0.25
K = GP.Mat32()
## Initialize model for the nonstationary properties GP
Random.seed!(13)
opt = TransD_GP.OptionsNonstat(optlog10λ,
                        nmin = nmin,
                        nmax = nmax,
                        fbounds = fbounds,
                        δ = δ,
                        demean = demean_ns,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        pnorm = pnorm,
                        K = K
                        )
calc_simple_RMS(img)

## set up McMC
nsamples, nchains, nchainsatone = 10001, 8, 1
Tmax = 2.50
addprocs(nchains)
@info "workers are $(workers())"
@everywhere any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
@everywhere using Distributed
@everywhere import MCMC_Driver
## run McMC
@time MCMC_Driver.main(optlog10λ, opt, img, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
#rmprocs(workers())
## plot
# MCMC_Driver.getchi2forall(opt_in)
# plot_last_target_model(img, opt_in)
