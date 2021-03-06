## stationary GP options
xall = transD_GP.get_image_prediction_points(img)
nmin, nmax= 2, 100
pnorm = 2.
K = transD_GP.GP.Mat32()
ftrain = img.d[img.select]
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
opt = transD_GP.OptionsStat(nmin = nmin,
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
