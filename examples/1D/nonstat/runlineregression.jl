using PyPlot, Random, Revise, Statistics, LinearAlgebra,
      Distributed, DelimitedFiles
srcdir = dirname(dirname(dirname(pwd())))*"/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
using GP, TransD_GP, GeophysOperator, MCMC_Driver
##1D functions
Random.seed!(10)
x = readdlm("func2.txt", ',', Float64, '\n')[:,1]
y = readdlm("func2.txt", ',', Float64, '\n')[:,2]
σ = 0.55
log10bounds = [-1.5 -0.5]
dec = 12
ynoisy = NaN .+ x
till = round(Int, length(x)/2)
ynoisy[1:till] = (σ*randn(size(y)) + y)[1:till]
ynoisy[till+1:dec:end] = (σ*randn(size(y)) + y)[till+1:dec:end]
line = GeophysOperator.Line(ynoisy;useML=false, σ=σ)
figure()
plot(x[:], y)
plot(x[:], ynoisy, ".m")
## make options for the multichannel lengthscale GP
δlog10λ = 0.2
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
fbounds = permutedims([extrema(ynoisy[.!isnan.(ynoisy)])...])
sdev_prop = 0.05*diff(fbounds, dims=2)[:]
sdev_pos = [0.05abs(diff([extrema(x)...])[1])]
demean_ns = true
K = GP.Mat32()
δ = 0.25
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
nsamples, nchains, nchainsatone = 50001, 4, 1
Tmax = 2.50
addprocs(nchains)
@info "workers are $(workers())"
@everywhere any($srcdir .== LOAD_PATH) || push!(LOAD_PATH, $srcdir)
@everywhere any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
@everywhere using Distributed
@everywhere import MCMC_Driver
## run McMC
@time MCMC_Driver.main(optlog10λ, opt, line, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
rmprocs(workers())
## plot
GeophysOperator.getchi2forall(opt)
ax = gcf().axes;
ax[3].set_ylim(100, 200)
ax[4].set_ylim(100, 200)
savefig("line_conv_ns_1.png", dpi=300)
GeophysOperator.getchi2forall(optlog10λ)
ax = gcf().axes;
ax[3].set_ylim(100, 200)
ax[4].set_ylim(100, 200)
savefig("line_conv_ns_2.png", dpi=300)
GeophysOperator.plot_posterior(line, opt, optlog10λ, burninfrac=0.2, figsize=(7.8,6), fsize=12)
ax = gcf().axes
ax[1].plot(ynoisy, x, ".c", alpha=0.4)
ax[1].plot(y, x, color="orange", alpha=0.4)
ax[1].plot(y, x, "--w", alpha=0.4)
ax[1].set_xlim(fbounds...)
savefig("jump1D_ns.png", dpi=300)
