#tests for "nuisance" inversions with HiQGA.transD_GP

## make options for a purely stationary properties GP
nmin, nmax = 2, 200
pnorm = 2.
λx, λy = 1, 1
x = 0:(0.005λx):λx
y = 0:(0.005λy):0.5λy
nmin, nmax = 2, 200
fbounds = [-2. 2]
δ = 0.1
sdev_prop = [0.1]
sdev_pos = [0.05, 0.05]
K = transD_GP.GP.Mat32()
λ = [0.05maximum(x), 0.05maximum(y)]
demean = false
sampledc = true
##
xall = zeros(2,length(x)*length(y))
for i in 1:size(xall,2)
    xid, yid = Tuple(CartesianIndices((length(x),length(y)))[i])
    xall[:,i] = [x[xid]; y[yid]]
end
xbounds = zeros(Float64,size(xall, 1), 2)
for dim in 1:size(xall, 1)
    xbounds[dim,:] = [minimum(xall[dim,:]), maximum(xall[dim,:])]
end

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
                        pnorm = pnorm,
                        needλ²fromlog = false,
                        updatenonstat = false,
                        quasimultid = false,
                        K = K,
                        timesλ = 3,
                        )

@time m = transD_GP.init(opt)

@testset "Creation tests for OptionsNuisance" begin
    nu_sdev = [0.1, 0.2, 0.3]
    nu_bounds = [0. 4.; -1.5 1.0; -1. 0.]
    optn = transD_GP.OptionsNuisance(opt;
        sdev = nu_sdev,
        bounds = nu_bounds,
        updatenuisances = true,
    )
    @testset "Test parameters set properly" begin
        @test isapprox(optn.sdev, [0.4, 0.5, 0.3], rtol=1e-4) #sdev in the struct is in "absolute units"
        @test all(optn.bounds .== nu_bounds)
        @test isapprox(optn.rotatebounds, nu_bounds .- mean(nu_bounds, dims=2), rtol=1e-4) #no covariance matrix was supplied
        @test optn.nnu == 3
        @test optn.updatenuisances == true
    end
    @testset "Constructor does not modify arguments" begin
        @test all(nu_sdev .== [0.1, 0.2, 0.3])
        @test all(nu_bounds .== [0. 4.; -1.5 1.0; -1. 0.])
    end
end