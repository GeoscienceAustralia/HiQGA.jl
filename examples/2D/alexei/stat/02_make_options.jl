## make the surface struct
using HiQGA.transD_GP, Random
surface = transD_GP.Surface(imagegrid;useML=false, σ=noisegrid)
## stationary GP options
nmin, nmax= 2, 100
pnorm = 2.
K = transD_GP.GP.Mat32()
ftrain = surface.d[surface.select]
fbounds = permutedims([extrema(ftrain)...])
δ = 0.3
λ = [50_000,50_000.]
demean = false
sampledc = true
sdev_prop = 0.05*diff(fbounds, dims=2)[:]
sdev_pos = 0.05*abs.([diff([extrema(xall[1,:])...])[1], diff([extrema(xall[2,:])...])[1]])
updatenonstat = false
needλ²fromlog = false
## Initialize a stationary GP using these options
Random.seed!(12)
opt = transD_GP.OptionsStat(nmin = nmin,
                        nmax = nmax,
                        xbounds = [minimum(xall[1,:]) maximum(xall[1,:]);
                                   minimum(xall[2,:]) maximum(xall[2,:])], 
                        fbounds = fbounds,
                        xall = xall,
                        λ = λ,
                        δ = δ,
                        demean = demean,
                        sampledc = sampledc,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        pnorm = pnorm,
                        quasimultid = false,
                        K = K,
                        timesλ = 3,
                        needλ²fromlog = needλ²fromlog,
                        updatenonstat = updatenonstat
                        )
