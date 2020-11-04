## make options for the multichannel lengthscale GP
log10bounds = [-1.7 -0.3]
δlog10λ = 0.2
nminlog10λ, nmaxlog10λ = 2, 30
pnorm = 2.
Klog10λ = transD_GP.GP.Mat32()
λlog10λ = [0.02abs(diff([extrema(x)...])[1])]
demean = false
sdev_poslog10λ = [0.05abs(diff([extrema(x)...])[1])]
sdev_proplog10λ = 0.05*diff(log10bounds, dims=2)[:]
xall = permutedims(collect(x))
xbounds = permutedims([extrema(x)...])
## Initialize a lengthscale model using these options
Random.seed!(12)
optlog10λ = transD_GP.OptionsStat(nmin = nminlog10λ,
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
                        save_freq = 20,
                        timesλ = 3.6,
                        peskycholesky = true
                        )
## make options for the nonstationary actual properties GP
nmin, nmax = 2, 30
fbounds = permutedims([extrema(ynoisy[.!isnan.(ynoisy)])...])
ymin, ymax = extrema(y)
fbounds[1] > ymin && (fbounds[1] = ymin)
fbounds[2] < ymax && (fbounds[2] = ymax)
sdev_prop = 0.05*diff(fbounds, dims=2)[:]
sdev_pos = [0.05abs(diff([extrema(x)...])[1])]
demean_ns = true
K = transD_GP.GP.Mat32()
δ = 0.2
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
nsamples, nchains, nchainsatone = 400001, 4, 1
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
    burninfrac=0.25, figsize=(7.5,4), fsize=8, nbins=100)
ax=gcf().axes
p = ax[1].scatter(ynoisy, x, c="w", alpha=0.2, s=25)
ax[1].plot(y, x, "--w", alpha=0.5)
del = fbounds[2]-fbounds[1]
ax[1].set_xlim(fbounds[1]-0.05del, fbounds[2]+0.05del,)
savefig("jump1D_ns.png", dpi=300)
ax[1].set_ylim(0.35,0.45)
p.set_alpha([0.5])
p.set_sizes([35])
ax[1].invert_yaxis()
savefig("jump1D_ns_zoom.png", dpi=300)
## correlated residuals (hacky, uncomment after running the stationary script first)
# h, edges = GeophysOperator.makehist(line, opt)
# ax[2].pcolormesh(x[sort(linidx)], edges, h)
# GeophysOperator.nicenup(f, fsize=8)
