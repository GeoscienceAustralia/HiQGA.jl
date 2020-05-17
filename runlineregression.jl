using PyPlot, Test, Random, Revise, Statistics, LinearAlgebra, Distributed
any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
using GP, TransD_GP, GeophysOperator, MCMC_Driver
## make options for the multichannel lengthscale GP
nminlog10λ, nmaxlog10λ = 2, 20
pnorm = 2.
Klog10λ = GP.Mat32()
λx = 1
x = 1:(0.01λx):2λx
λlog10λ = [0.05abs(diff([extrema(x)...])[1])]
demean = false
sdev_poslog10λ = [0.01maximum(x)]
log10bounds = [-1 -0.69]
δlog10λ = 0.1
sdev_proplog10λ = [0.1]
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
fbounds = [-2. 2]
δ = 0.5
sdev_prop = [0.1]
sdev_pos = [0.05]
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
##
nsamples, nchains, nchainsatone = 500001, 4, 1
Tmax = 2.50
line = Line([1])
addprocs(nchains)
@info "workers are $(workers())"
@everywhere any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
@everywhere using Distributed
@everywhere import MCMC_Driver
## run McMC
@time MCMC_Driver.main(optlog10λ, opt, line, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
rmprocs(workers())
