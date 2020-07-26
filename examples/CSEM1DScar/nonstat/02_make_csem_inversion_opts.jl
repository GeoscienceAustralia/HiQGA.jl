using GP, TransD_GP, GeophysOperator, MCMC_Driver, Distributed
## make options for the multichannel lengthscale GP
nminlog10λ, nmaxlog10λ = 2, 40
pnorm = 2.
Klog10λ = GP.Mat32()
log10bounds = [0 1.5]
λlog10λ = [0.07abs(diff([extrema(znall)...])[1])]
δlog10λ = 0.1
demean = false
sdev_poslog10λ = [0.035abs(diff([extrema(znall)...])[1])]
sdev_proplog10λ = 0.05*diff(log10bounds, dims=2)[:]
xall = permutedims(collect(znall))
xbounds = permutedims([extrema(znall)...])
## Initialize a lengthscale model using these options
Random.seed!(12)
optlog10λ = TransD_GP.OptionsStat(nmin = nminlog10λ,
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
                        K = Klog10λ,
                        timesλ = 3.6
                        )
## make options for the nonstationary actual properties GP
nmin, nmax = 2, 20
fbounds = [-0.5 2.3]
δ = 0.2
sdev_prop = 0.01*diff(fbounds, dims=2)[:]
sdev_pos = [0.006abs(diff([extrema(znall)...])[1])]
demean_ns = true
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
