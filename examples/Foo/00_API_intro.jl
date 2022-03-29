using HiQGA.transD_GP

# Forward setttings
# dummy physics struct defined in src/FooPhysics.jl
F = transD_GP.FooPhysics.Foo(25, rand(10))

# some dummy forward
@info transD_GP.FooPhysics.returnphysics!(F, rand(3))

# struct for data and physics operator defined in src/FooPhysicsInversion.jl 
FOp = transD_GP.FooPhysicsInversion.FooInversion(rand(10), F)

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
