using .transD_GP, Distributed
## make options for the stationary GP
fileprefix = "OrstUhn_"
K = transD_GP.GP.OrstUhn()
nmin, nmax = 2, 50
fbounds = [-0.5 3.25]
sdev_pos = [0.05abs(diff([extrema(znall)...])[1])]
sdev_prop = 0.07*diff(fbounds, dims=2)[:]
demean = false
sdev_dc = 0.008*diff(fbounds, dims=2)[:]
sampledc = true
xall = permutedims(collect(znall))
xbounds = permutedims([extrema(znall)...])
λ, δ = [2], 0.1
F.Fvec[1].useML = false
F.Fvec[2].useML = false

opt = transD_GP.OptionsStat(nmin = nmin,
                        nmax = nmax,
                        xbounds = xbounds,
                        fbounds = fbounds,
                        fdataname = fileprefix,
                        xall = xall,
                        λ = λ,
                        δ = δ,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        save_freq = 25,
                        demean = demean,
                        sampledc = sampledc,
                        sdev_dc = sdev_dc,
                        quasimultid = false,
                        K = K
                        )
