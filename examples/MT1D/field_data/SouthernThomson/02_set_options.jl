## make transD options
# prior stuff
nmin, nmax = 2, 50
K = transD_GP.GP.OrstUhn()
fbounds = [-1 3.0]
xall = permutedims(collect(znall)) # discretization 
xbounds = permutedims([extrema(znall)...]) # discretization bounds
xall = permutedims(collect(znall))
xbounds = permutedims([extrema(znall)...])
λ, δ = [2], 0.1
# proposal stuff
sdev_pos = [0.05abs(diff([extrema(znall)...])[1])]
sdev_prop = 0.07*diff(fbounds, dims=2)[:]
sdev_dc = 0.01*diff(fbounds, dims=2)[:]
## initialise options
opt = transD_GP.OptionsStat(;
                        nmin = nmin,
                        nmax = nmax,
                        xbounds = xbounds,
                        fbounds = fbounds,
                        fdataname = fname*"_",
                        xall = xall,
                        λ = λ,
                        δ = δ,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        sdev_dc = sdev_dc,
                        save_freq = 25,
                        quasimultid = false,
                        K = K
                        )