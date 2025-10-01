## make options for the stationary GP
fileprefix = "SWD_test_"
K = transD_GP.GP.OrstUhn()
nmin, nmax = 2, 40
fbounds = [minimum(vs) maximum(vs)]
sdev_pos = [0.05abs(diff([extrema(znall)...])[1])]
sdev_prop = 0.07*diff(fbounds, dims=2)[:]
demean = false
sdev_dc = 0.008*diff(fbounds, dims=2)[:]
sampledc = true
xall = permutedims(collect(znall))
xbounds = permutedims([extrema(znall)...])
λ, δ = [5], 0.1
## Initialize a stationary GP using these options
using Random
Random.seed!(12)
opt = transD_GP.OptionsStat(;nmin, nmax, xbounds, fbounds,
                        fdataname = fileprefix, xall, λ, δ, sdev_prop, sdev_pos,
                        demean, sampledc, sdev_dc, K,
                        quasimultid = false, save_freq = 50, fdataname = fileprefix,
                        )
## 
transD_GP.get_misfit(vs, opt, swd)
