## make options for the stationary GP
fileprefix = "SWD_test_"
K = transD_GP.GP.OrstUhn() # GP Kernel type, OrstUhn is L1 sharp
nmin, nmax = 2, 40 # no. of GP nuclei
fbounds = [0.9minimum(vs) 1.1maximum(vs)] # property bounds, n_property × 2 array
xall = permutedims(collect(znall)) # spatial locations to interpolate GP at
xbounds = permutedims([extrema(znall)...]) # spatial bounds, n_spatial_dims × 2 array
sdev_pos = [0.05abs(diff([extrema(znall)...])[1])] # perturbation for position, n_spatial_dims × 2 array
sdev_prop = 0.07*diff(fbounds, dims=2)[:] # perturbation for property, n_property × 2 array
sdev_dc = 0.008*diff(fbounds, dims=2)[:] # DC perturbation for property, n_property × 2 array
λ, δ = [5.5], 0.1 # correlation length, nugget
## Initialize a stationary GP using these options
using Random
Random.seed!(12)
opt = transD_GP.OptionsStat(;nmin, nmax, xbounds, fbounds,
                        xall, λ, δ, sdev_prop, sdev_pos,
                        demean, sampledc, sdev_dc, K,
                        quasimultid = false, 
                        save_freq = 25, 
                        stat_window = 200,
                        fdataname = fileprefix,
                        )
## 
@info "@info χ² is $(transD_GP.get_misfit(vs, opt, swd))" # check χ²
