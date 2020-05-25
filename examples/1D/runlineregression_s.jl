using PyPlot, Test, Random, Revise, Statistics, LinearAlgebra, Distributed
any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
using GP, TransD_GP, GeophysOperator, MCMC_Driver, DelimitedFiles
##1D functions
easy = false
Random.seed!(10)
fractrain = 1
if easy
    Random.seed!(10)
    x = LinRange(-2,2,101)
    y = sin.(x) +2exp.(-30x.^2)
    σ = 0.3
    δ = 0.25
else
    Random.seed!(10)
    x = LinRange(0, 1, 201)
    y = readdlm("func.txt")[:]
    σ = 0.55
    δ = 0.25
    λ = [0.05]
end
ntrain = round(Int, (1-fractrain)*length(y))
linidx = randperm(length(y))[1:ntrain]
ynoisy = σ*randn(size(y)) + y
ynoisy[linidx] .= NaN
line = GeophysOperator.Line(ynoisy;useML=false, σ=σ)
figure()
plot(x[:], y)
plot(x[:], ynoisy, ".m")
## make options for the stationary GP
nmin, nmax = 2, 20
pnorm = 2.
K = GP.Mat32()
demean = true
fbounds = permutedims([extrema(ynoisy[.!isnan.(ynoisy)])...])
sdev_pos = [0.05abs(diff([extrema(x)...])[1])]
sdev_prop = 0.05*diff(fbounds, dims=2)[:]
xall = permutedims(collect(x))
xbounds = permutedims([extrema(x)...])
updatenonstat = false
needλ²fromlog = false
## Initialize a lengthscale model using these options
Random.seed!(12)
opt = TransD_GP.OptionsStat(nmin = nmin,
                        nmax = nmax,
                        xbounds = xbounds,
                        fbounds = fbounds,
                        xall = xall,
                        λ = λ,
                        δ = δ,
                        demean = demean,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        pnorm = pnorm,
                        quasimultid = false,
                        K = K,
                        timesλ = 3.,
                        needλ²fromlog = needλ²fromlog,
                        updatenonstat = updatenonstat
                        )
## Initialize options for the dummy nonstationary properties GP
Random.seed!(13)
optdummy = TransD_GP.OptionsNonstat(opt,
                        nmin = nmin,
                        nmax = nmax,
                        fbounds = fbounds,
                        δ = δ,
                        demean = demean,
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
@time MCMC_Driver.main(opt, optdummy, line, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
rmprocs(workers())
## plot
GeophysOperator.getchi2forall(opt)
GeophysOperator.plot_posterior(line, opt, burninfrac=0.2, figsize=(5,10))
ax = gcf().axes
ax[1].plot(ynoisy, x, ".m")
ax[1].plot(y, x, color="orange")
ax[1].plot(y, x, "--k")
