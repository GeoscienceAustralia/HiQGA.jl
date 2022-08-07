using .transD_GP, Distributed
## make options for the stationary GP
fileprefix = "transD_"
# prior specifications
K = transD_GP.GP.OrstUhn() # type of GP kernel
nmin, nmax = 2, 50 # number of GP nuclei
fbounds = [-0.5 3.25] # log10 RESISTIVITY bounds to sample between
xall = permutedims(znall) # depth locations to interpolate resistivity to 
xbounds = permutedims([extrema(znall)...]) # bounds of the depth locations
λ, δ = [2], 0.1 # correlation length in depth units, tolerance for GP resistivity
# proposal specifications computed as fractions of prior bounds
sdev_pos = 0.05vec(diff([extrema(znall)...]))
sdev_prop = [0.07*diff(fbounds, dims=2)...]
sdev_dc = [0.008*diff(fbounds, dims=2)...]
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
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        sdev_dc = sdev_dc,
                        quasimultid = false,
                        save_freq = 50,
                        K = K
                        )
