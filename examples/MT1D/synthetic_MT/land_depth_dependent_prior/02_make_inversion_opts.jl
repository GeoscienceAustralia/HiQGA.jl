## make options for the stationary GP
fileprefix = "test_stretchis$(string(F.stretch))_"
# McMC priors
K = transD_GP.GP.OrstUhn()
nmin, nmax = 2, 40 # no. of nuclei
fbounds = F.stretch ? [0 1.] : [-0.5 3] # stretch priors or MT bounds in log10rho
xall = permutedims(collect(znall)) # discretization 
xbounds = permutedims([extrema(znall)...]) # discretization bounds
λ, δ = [2], F.stretch ? 0.03 : 0.1 # correlation length units, nugget in stretch units or log10 resistivity
# McMC proposals
sdev_pos = [0.05abs(diff([extrema(znall)...])[1])]
sdev_prop = 0.07*diff(fbounds, dims=2)[:]
sdev_dc = 0.008*diff(fbounds, dims=2)[:]
## Initialize a stationary GP using these options
using Random
Random.seed!(12)
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
