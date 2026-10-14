transD_GP.getchi2forall(opt)
##
opt.xall[:] .= zall
transD_GP.plot_posterior(tempest, opt, burninfrac=0.5, figsize=(6,6), cmappdf="hot",
        qp1=0.05, qp2=0.95, nbins=200, vmaxpc=1)
ax = gcf().axes
ax[1].invert_xaxis()
ax[1].step(log10.(ρ[2:end]), z[2:end], color="r", alpha=0.5)
ax[1].step(log10.(ρ[2:end]), z[2:end], color="m", alpha=0.5, linestyle="--")
ax[1].set_ylim(280,0)
## nuisance histograms
transD_GP.plot_posterior(tempest, optn, burninfrac=0.5, nbins=50)
