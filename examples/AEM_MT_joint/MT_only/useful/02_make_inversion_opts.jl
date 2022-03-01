## make options for the stationary GP
fileprefix = "MT_only_useful_stretchis$(string(F.stretch))_"
K = transD_GP.GP.OrstUhn()
nmin, nmax = 2, 40
fbounds = F.stretch ? [0 1.] : [-0.5 4]
sdev_pos = [0.05abs(diff([extrema(znall)...])[1])]
sdev_prop = 0.07*diff(fbounds, dims=2)[:]
demean = false
sdev_dc = 0.008*diff(fbounds, dims=2)[:]
sampledc = true
xall = permutedims(collect(znall))
xbounds = permutedims([extrema(znall)...])
λ, δ = [2], F.stretch ? 0.03 : 0.1 
debug = false
## Initialize a stationary GP using these options
using Random
Random.seed!(12)
opt = transD_GP.OptionsStat(debug = debug,
                        nmin = nmin,
                        nmax = nmax,
                        xbounds = xbounds,
                        fbounds = fbounds,
                        fdataname = fileprefix,
                        xall = xall,
                        λ = λ,
                        δ = δ,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        save_freq = 50,
                        demean = demean,
                        sampledc = sampledc,
                        sdev_dc = sdev_dc,
                        quasimultid = false,
                        K = K
                        )
