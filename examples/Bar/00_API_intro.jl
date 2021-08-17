using transD_GP
include("BarPhysics.jl")
include("BarPhysicsInversion.jl")
using .BarPhysicsInversion
# Forward setttings
# dummy physics struct defined in BarPhysics.jl
F = BarPhysics.Bar(25, rand(10))

# some dummy forward
@info BarPhysics.returnphysics!(F, rand(3))

# struct for data and physics operator defined in BarPhysicsInversion.jl 
FOp = BarPhysicsInversion.BarInversion(rand(10), F)

## Inverse settings
# prior bounds and setings
fbounds = [0 1.] 
xall = collect(permutedims(0:0.1:1))
xbounds = permutedims([extrema(xall)...])
λ = vec(diff(xbounds, dims=2))
nmin, nmax = 2, 5
demean, sampledc = false, true
# McMC proposals
δ = 0.1*diff(fbounds, dims=2)[1]
sdev_pos = abs.(vec(diff(xbounds, dims=2)))
sdev_prop = abs.(vec(0.05*diff(fbounds, dims=2)))
K = transD_GP.GP.OrstUhn()

# make MCMC options using prior and proposals
opt = transD_GP.OptionsStat(nmin = nmin,
                        nmax = nmax,
                        xbounds = xbounds,
                        fbounds = fbounds,
                        xall = xall,
                        λ = λ,
                        δ = δ,
                        demean = demean,
                        sampledc = sampledc,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        quasimultid = false,
                        K = K
                        )

## Put forward and inverse together
# initialise a random model
m = transD_GP.init(opt) 

# compute misfit
transD_GP.get_misfit(m, opt, FOp)

# if you've gotten till here you're golden
# you don't really need to initialize a model and compute
# misfit outside a module