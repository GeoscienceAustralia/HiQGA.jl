using transD_GP, Distributed
## make options for the stationary GP
fileprefix = "OrstUhn_"
K = transD_GP.GP.OrstUhn()
nmin, nmax = 2, 40
fbounds = [-0.5 2.5]
sdev_pos = [0.05abs(diff([extrema(znall)...])[1])]
sdev_prop = 0.05*diff(fbounds, dims=2)[:]
xall = permutedims(collect(znall))
xbounds = permutedims([extrema(znall)...])
λ, δ = [2], 0.05
## Initialize a stationary GP using these options
Random.seed!(12)
opt = transD_GP.OptionsStat(nmin = nmin,
                        nmax = nmax,
                        xbounds = xbounds,
                        fbounds = fbounds,
                        fdataname = fileprefix,
                        xall = xall,
                        λ = λ,
                        δ = δ,
                        demean = demean,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        save_freq = 25,
                        quasimultid = false,
                        K = K
                        )
