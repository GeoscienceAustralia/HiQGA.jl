## make options for the stationary GP
λ = [50.]
δ = 0.2
nmin, nmax = 2, 100
pnorm = 2.
K = transD_GP.GP.Mat32()
demean = true
fbounds = permutedims([extrema(ynoisy[.!isnan.(ynoisy)])...])
ymin, ymax = extrema(y)
fbounds[1] > ymin && (fbounds[1] = ymin)
fbounds[2] < ymax && (fbounds[2] = ymax)
sdev_pos = [0.05abs(diff([extrema(x)...])[1])]
sdev_prop = 0.05*diff(fbounds, dims=2)[:]
xall = permutedims(collect(x))
xbounds = permutedims([extrema(x)...])
updatenonstat = false
needλ²fromlog = false
## Initialize a stationary GP using these options
Random.seed!(12)
opt = transD_GP.OptionsStat(nmin = nmin,
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
                        timesλ = 3.6,
                        save_freq = 20,
                        needλ²fromlog = needλ²fromlog,
                        updatenonstat = updatenonstat,
                        peskycholesky = true,
                        )
## set up McMC
nsamples, nchains, nchainsatone = 10001, 8, 1
Tmax = 2.50
addprocs(nchains)
@info "workers are $(workers())"
@everywhere using Distributed
@everywhere using transD_GP
## run McMC
@time transD_GP.main(opt, line, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
rmprocs(workers())
## plot
transD_GP.getchi2forall(opt, fsize=8, alpha=0.5)
transD_GP.plot_posterior(line, opt, burninfrac=0.5, figsize=(8.5,4), fsize=8, nbins=100)
ax = gcf().axes
p = ax[1].scatter(ynoisy, x, c="m", alpha=0.2, s=2)
# ax[1].plot(y, x, "--w", alpha=0.5)
del = fbounds[2]-fbounds[1]
ax[1].set_xlim(fbounds[1]-0.05del, fbounds[2]+0.05del,)
savefig("jump1D_high.png", dpi=300)
## zoom in
ax[1].set_ylim(0.35,0.45)
p.set_alpha(0.5)
p.set_sizes([35])
gcf().text(0.02, 0.9, "a.", fontsize=14, color="red")
ax[1].invert_yaxis()
savefig("jump1D_high_zoom.png", dpi=300)
## residuals jump at 0.4
# h, edges = transD_GP.makehist(line, opt)
# f, ax = plt.subplots(2,1, figsize=(4,3), sharex=true, sharey=true)
# ax[1].pcolormesh(x[sort(linidx)], edges, h)
# ax[1].set_ylabel("residual")
# ax[2].set_ylabel("residual")
# ax[1].set_title("Stationary TDGP")
# ax[2].set_xlabel("depth m")
# ax[2].set_title("Non-stationary TDGP")
# gcf()
