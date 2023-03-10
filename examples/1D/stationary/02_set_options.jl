## make options for the stationary GP
# prior stuff
λ = [.1] # correlation length
δ = 0.05 # nugget
nmin, nmax = 2, 30 # no. of nuclei
K = transD_GP.GP.Mat32() # Kernel function 
fbounds = permutedims([extrema(ynoisy[.!isnan.(ynoisy)])...]) # GP value bounds
ymin, ymax = extrema(y)
fbounds[1] > ymin && (fbounds[1] = ymin)
fbounds[2] < ymax && (fbounds[2] = ymax)
xall = permutedims(collect(x)) # discretization 
xbounds = permutedims([extrema(x)...]) # discretization bounds
# McMC proposal stuff
sdev_pos = [0.05abs(diff([extrema(x)...])[1])]
sdev_prop = 0.05*diff(fbounds, dims=2)[:]
## Initialize a stationary GP using these options
opt = transD_GP.OptionsStat(nmin = nmin,
                        nmax = nmax,
                        xbounds = xbounds,
                        fbounds = fbounds,
                        xall = xall,
                        λ = λ,
                        δ = δ,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        quasimultid = false,
                        K = K,
                        save_freq = 20,
                        )
