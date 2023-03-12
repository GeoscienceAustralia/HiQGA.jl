using .transD_GP, Distributed
## make options for the stationary GP
fileprefix = "JointInv_"
# McMC priors
K = transD_GP.GP.OrstUhn()
nmin, nmax = 2, 50 # no. of nuclei
fbounds = [-0.5 3.25] # bounds in log 10 rho
xall = permutedims(collect(znall)) # discretization
xbounds = permutedims([extrema(znall)...]) # discretization bounds
λ, δ = [2], 0.1
# McMC proposal stuff
sdev_pos = [0.05abs(diff([extrema(znall)...])[1])]
sdev_prop = 0.07*diff(fbounds, dims=2)[:]
sdev_dc = 0.008*diff(fbounds, dims=2)[:]
# for maximum likelihood errors
F.Fvec[1].useML = false
F.Fvec[2].useML = false
## Make a stationary GP with above options
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
