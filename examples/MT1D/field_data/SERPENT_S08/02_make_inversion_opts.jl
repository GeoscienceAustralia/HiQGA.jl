## make options for the stationary GP
# prior stuff
fileprefix = "MT_SERPENT_08_stretchis$(string(F.stretch))_"
K = transD_GP.GP.OrstUhn()
nmin, nmax = 2, 40
fbounds = F.stretch ? [0 1.] : [-1 5]
xall = permutedims(collect(znall))
xbounds = permutedims([extrema(znall)...])
λ, δ = [2], F.stretch ? 0.02 : 0.1
# proposal stuff
sdev_pos = [0.05abs(diff([extrema(znall)...])[1])]
sdev_prop = 0.07*diff(fbounds, dims=2)[:]
sdev_dc = 0.01*diff(fbounds, dims=2)[:]
## Initialize a stationary GP using these options
opt = transD_GP.OptionsStat(;
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
                        save_freq = 25,
                        sdev_dc = sdev_dc,
                        quasimultid = false,
                        K = K
                        )
