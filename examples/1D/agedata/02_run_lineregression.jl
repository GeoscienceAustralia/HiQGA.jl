## make options for the lengthscale GP
log10bounds = [1.3 3]
δlog10λ = 0.075
nminlog10λ, nmaxlog10λ = 2, 300
pnorm = 2.
Klog10λ = transD_GP.GP.OrstUhn()
λlog10λ = [0.03abs(diff([extrema(X)...])[1])]
demean = false
sdev_poslog10λ = [0.05abs(diff([extrema(X)...])[1])]
sdev_proplog10λ = 0.1*diff(log10bounds, dims=2)[:]
xall = permutedims(collect(X))
xbounds = permutedims([extrema(X)...])
## Initialize a lengthscale model using these options
Random.seed!(12)
optlog10λ = transD_GP.OptionsStat(nmin = nminlog10λ,
                        nmax = nmaxlog10λ,
                        xbounds = xbounds,
                        fbounds = log10bounds,
                        fdataname = "variable_",
                        xall = xall,
                        λ = λlog10λ,
                        δ = δlog10λ,
                        demean = demean,
                        sdev_prop = sdev_proplog10λ,
                        sdev_pos = sdev_poslog10λ,
                        pnorm = pnorm,
                        quasimultid = false,
                        K = Klog10λ,
                        save_freq = 100,
                        timesλ = 4,
                        peskycholesky = true
                        )
## make options for the nonstationary actual properties GP
nmin, nmax = 2, 100
fbounds = permutedims([extrema(Y[.!isnan.(Y)])...])
ymin, ymax = extrema(Y)
fbounds[1] > ymin && (fbounds[1] = ymin)
fbounds[2] < ymax && (fbounds[2] = ymax)
sdev_prop = 0.015*diff(fbounds, dims=2)[:]
sdev_pos = [0.01abs(diff([extrema(X)...])[1])]
demean_ns = true
K = transD_GP.GP.OrstUhn()
δ = 0.075
## Initialize model for the nonstationary properties GP
Random.seed!(13)
opt = transD_GP.OptionsNonstat(optlog10λ,
                        nmin = nmin,
                        nmax = nmax,
                        fbounds = fbounds,
                        δ = δ,
                        demean = demean_ns,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        pnorm = pnorm,
                        K = K,
                        )
## set up McMC
nsamples, nchains, nchainsatone = 400001, 8, 1
Tmax = 2.50
addprocs(nchains)
@info "workers are $(workers())"
@everywhere using Distributed
@everywhere using transD_GP
## run McMC
@time transD_GP.main(optlog10λ, opt, line, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
rmprocs(workers())
## plot
transD_GP.getchi2forall(opt, fsize=8, alpha=0.5)
transD_GP.getchi2forall(optlog10λ, fsize=8, alpha=0.5)
transD_GP.plot_posterior(line, opt, optlog10λ,
    burninfrac=0.5, figsize=(11.5,11.5), fsize=8, nbins=100, vmaxpc=0.5, cmappdf="bone")
ax=gcf().axes
p = ax[1].scatter(Y, X, c="orange", alpha=0.05, s=2)
del = fbounds[2]-fbounds[1]
ax[1].set_xlim(fbounds[1]-0.05del, fbounds[2]+0.05del,)
savefig("jump1D_ns.png", dpi=300)
