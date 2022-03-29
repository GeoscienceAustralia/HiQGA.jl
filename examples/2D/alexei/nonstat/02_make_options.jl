## make the surface struct
using HiQGA.transD_GP, Random
surface = transD_GP.Surface(imagegrid;useML=false, σ=noisegrid)
## lengthscale GP options (stationary)
nminlog10λ, nmaxlog10λ = 2, 150
pnorm = 2.
Klog10λ = transD_GP.GP.Mat32()
log10bounds = [log10(20_000) log10(300_000); log10(20_000) log10(300_000)]
δlog10λ = 0.2
# making sure it is not anisotropic ...
λlog10λ = 0.1*abs.([diff([extrema(xall[2,:])...])[1], diff([extrema(xall[2,:])...])[1]])
demean = false # essential for length scale change approximations
sdev_proplog10λ = 0.05*diff(log10bounds, dims=2)[:]
sdev_poslog10λ = 0.05*abs.([diff([extrema(xall[1,:])...])[1], diff([extrema(xall[2,:])...])[1]])
## Initialize a lengthscale model using these options
Random.seed!(12)
optlog10λ = transD_GP.OptionsStat(nmin = nminlog10λ,
                        nmax = nmaxlog10λ,
                        xbounds = [minimum(xall[1,:]) maximum(xall[1,:]);
                                   minimum(xall[2,:]) maximum(xall[2,:])], 
                        fbounds = log10bounds,
                        xall = xall,
                        λ = λlog10λ,
                        δ = δlog10λ,
                        needλ²fromlog = true,
                        updatenonstat = true,
                        demean = demean,
                        sdev_prop = sdev_proplog10λ,
                        sdev_pos = sdev_poslog10λ,
                        pnorm = pnorm,
                        quasimultid = false,
                        K = Klog10λ,
                        timesλ = 4.
                        )
## make options for the nonstationary actual properties GP
nmin, nmax = 2, 150
ftrain = surface.d[surface.select]
fbounds = permutedims([extrema(ftrain)...])
sdev_prop = 0.05*diff(fbounds, dims=2)[:]
sdev_pos = 0.05*abs.([diff([extrema(xall[1,:])...])[1], diff([extrema(xall[2,:])...])[1]])
demean_ns = false
sampledc = true
δ = 1.
K = transD_GP.GP.Mat32()
## Initialize model for the nonstationary properties GP
Random.seed!(13)
opt = transD_GP.OptionsNonstat(optlog10λ,
                        nmin = nmin,
                        nmax = nmax,
                        fbounds = fbounds,
                        δ = δ,
                        demean = demean_ns,
                        sampledc = sampledc,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        pnorm = pnorm,
                        K = K
                        )
