## lengthscale GP options
xall = transD_GP.get_image_prediction_points(img)
nminlog10λ, nmaxlog10λ = 2, 100
pnorm = 2.
Klog10λ = transD_GP.GP.Mat32()
log10bounds = [log10(10) log10(2000); log10(10) log10(2000)]
δlog10λ = 0.3
λlog10λ = 0.1*abs.([diff([extrema(xall[1,:])...])[1], diff([extrema(xall[2,:])...])[1]])
demean = false # essential for length scale change approximations
sdev_proplog10λ = 0.05*diff(log10bounds, dims=2)[:]
sdev_poslog10λ = 0.05*abs.([diff([extrema(xall[1,:])...])[1], diff([extrema(xall[2,:])...])[1]])
## Initialize a lengthscale model using these options
Random.seed!(12)
optlog10λ = transD_GP.OptionsStat(nmin = nminlog10λ,
                        nmax = nmaxlog10λ,
                        xbounds = [img.x[1] img.x[end];img.y[1] img.y[end]],
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
                        timesλ = 3.6
                        )
## make options for the nonstationary actual properties GP
nmin, nmax = 2, 100
ftrain = img.d[img.select]
fbounds = permutedims([extrema(ftrain)...])
sdev_prop = 0.05*diff(fbounds, dims=2)[:]
sdev_pos = 0.05*abs.([diff([extrema(xall[1,:])...])[1], diff([extrema(xall[2,:])...])[1]])
demean_ns = true
δ = 0.25
K = transD_GP.GP.Mat32()
## Initialize model for the nonstationary properties GP
Random.seed!(13)
opt = transD_GP.OptionsNonstat(optlog10λ,
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
