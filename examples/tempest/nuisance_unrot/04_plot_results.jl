transD_GP.getchi2forall(optn)
gcf()
##
opt.xall[:] .= zall
transD_GP.plot_posterior(tempest, opt, burninfrac=0.5, figsize=(10,6),
    cmappdf="inferno", qp1=0.05, qp2=0.95, nbins=200, vmaxpc=0.5)
ax = gcf().axes
ax[1].invert_xaxis()
ax[1].step(log10.(ρ[2:end]), z[2:end], color="w", alpha=0.8)
## nuisance histograms
transD_GP.plot_posterior(tempest, optn, burninfrac=0.5, nbins=50)
gcf()
