## make options for the stationary GP
fileprefix = "MT_synth_"
# McMC prior stuff
K = transD_GP.GP.OrstUhn()
nmin, nmax = 2, 40 # no. of l=nuclei
fbounds = [-1 3] # MT bounds in log 10 rho
xall = permutedims(collect(znall)) # discretization
xbounds = permutedims([extrema(znall)...]) # discretization bounds
λ, δ = [2], 0.2 # correlation length units, nugget in log10 resistivity
# McMC proposal stuff
sdev_pos = [0.05abs(diff([extrema(znall)...])[1])]
sdev_prop = 0.07*diff(fbounds, dims=2)[:]
sdev_dc = 0.008*diff(fbounds, dims=2)[:]
## Initialize a stationary GP using these options
using Random
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
                        sdev_dc = sdev_dc,
                        quasimultid = false,
                        K = K
                        )
