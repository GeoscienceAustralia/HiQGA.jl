## make options for the stationary GP
λ = [2]
δ = 0.1
nmin, nmax = 2, 40
K = transD_GP.GP.OrstUhn()
demean = false
sampledc = true
fbounds = [lo hi]
ρmin, ρmax = extrema(ρ)
sdev_pos = [0.05abs(diff([extrema(znall)...])[1])]
sdev_prop = 0.05*diff(fbounds, dims=2)[:]
xall = permutedims(collect(znall))
xbounds = permutedims([extrema(znall)...])
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
                        sampledc = sampledc,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        K = K,
                        save_freq = 50,
                        quasimultid = false
                        )
## set up McMC
using Distributed
nsamples, nchains, nchainsatone = 100_001, 4, 1
Tmax = 2.50
addprocs(nchains)
@info "workers are $(workers())"
@everywhere using Distributed 
@everywhere using HiQGA.transD_GP
## run McMC
@time transD_GP.main(opt, line, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
rmprocs(workers())
## plot
transD_GP.getchi2forall(opt, fsize=8, alpha=0.5)
##
opt.xall[:] = zall[:]
c1, cm, c3, =transD_GP.plot_posterior(line, opt, showslope=false, plotmean=false, qp1=0.1, qp2=0.9, 
    plotCI=true,pdfnormalize=false, burninfrac=0.1, figsize=(6,8), fsize=12, nbins=50, lwidth=0.5)
ax = gcf().axes
p = ax[1].scatter(log10.(ρnoisy), zall, c="w", alpha=0.5, s=50)
p = ax[1].scatter(log10.(ρnoisy), zall, c="g", alpha=0.5, s=30)
ax[1].step(log10.(ρ), z, "-w", linewidth=2)
ax[1].step(log10.(ρ), z, "--k", linewidth=1)
del = fbounds[2]-fbounds[1]
ax[1].set_xlim(fbounds[1]-0.05del, fbounds[2]+0.05del,)
ax[1].invert_xaxis()
ax = ax[1].axis()
##
f = figure(figsize=(6,8))
step(m_ridge, z, linewidth=1, label="ridge")
step(m_mle, z, linewidth=2, color="m", label="MLE")
step(ocm, z, linewidth=2, color="r", label="occam")
step(log10.(ρ), z, color="k", linewidth=2, label="true")
fill_betweenx(zall, c1, c3, alpha=0.25)
gca().axis(ax)
xlabel(L"\log_{10}\rho")
ylabel("depth m")
transD_GP.nicenup(gcf(), fsize=12)