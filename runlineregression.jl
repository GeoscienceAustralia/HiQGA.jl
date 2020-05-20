using PyPlot, Test, Random, Revise, Statistics, LinearAlgebra, Distributed
any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
using GP, TransD_GP, GeophysOperator, MCMC_Driver, DelimitedFiles
##1D functions
easy = false
if easy
    Random.seed!(10)
    x = LinRange(-2,2,101)
    y = sin.(x) +2exp.(-30x.^2)
    σ = 0.3
    log10bounds = [-0.5 0.1]
    δlog10λ = 0.3
    δ = 0.25
else
    Random.seed!(10)
    x = LinRange(0, 1, 201)
    y = readdlm("func.txt")[:]
    σ = 0.55
    log10bounds = [-1.5 -0.5]
    δlog10λ = 0.2
    δ = 0.25
end
ynoisy = σ*randn(size(y)) + y
line = GeophysOperator.Line(ynoisy;useML=false, σ=σ)
figure()
plot(x[:], y)
plot(x[:], ynoisy, ".m")
## make options for the multichannel lengthscale GP
nminlog10λ, nmaxlog10λ = 2, 20
pnorm = 2.
Klog10λ = GP.Mat32()
λlog10λ = [0.05abs(diff([extrema(x)...])[1])]
demean = false
sdev_poslog10λ = [0.05abs(diff([extrema(x)...])[1])]
sdev_proplog10λ = 0.05*diff(log10bounds, dims=2)[:]
xall = permutedims(collect(x))
xbounds = permutedims([extrema(x)...])
## Initialize a lengthscale model using these options
Random.seed!(12)
optlog10λ = TransD_GP.OptionsStat(nmin = nminlog10λ,
                        nmax = nmaxlog10λ,
                        xbounds = xbounds,
                        fbounds = log10bounds,
                        xall = xall,
                        λ = λlog10λ,
                        δ = δlog10λ,
                        demean = demean,
                        sdev_prop = sdev_proplog10λ,
                        sdev_pos = sdev_poslog10λ,
                        pnorm = pnorm,
                        quasimultid = false,
                        K = Klog10λ,
                        timesλ = 3.
                        )
## make options for the nonstationary actual properties GP
nmin, nmax = 2, 20
fbounds = permutedims([extrema(ynoisy)...])
sdev_prop = 0.05*diff(fbounds, dims=2)[:]
sdev_pos = [0.05abs(diff([extrema(x)...])[1])]
demean_ns = true
K = GP.Mat32()
## Initialize model for the nonstationary properties GP
Random.seed!(13)
opt = TransD_GP.OptionsNonstat(optlog10λ,
                        nmin = nmin,
                        nmax = nmax,
                        fbounds = fbounds,
                        δ = δ,
                        demean = demean_ns,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        pnorm = pnorm,
                        K = K
                        )
## set up McMC
nsamples, nchains, nchainsatone = 15001, 4, 1
Tmax = 2.50
addprocs(nchains)
@info "workers are $(workers())"
@everywhere any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
@everywhere using Distributed
@everywhere import MCMC_Driver
## run McMC
@time MCMC_Driver.main(optlog10λ, opt, line, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
rmprocs(workers())
## plot
GeophysOperator.getchi2forall(opt)
GeophysOperator.getchi2forall(optlog10λ)
GeophysOperator.plot_posterior(line, opt, optlog10λ, burninfrac=0.2)
ax = gcf().axes
ax[1].plot(ynoisy, x, ".m")
ax[1].plot(y, x, color="orange")
ax[1].plot(y, x, "--k")
