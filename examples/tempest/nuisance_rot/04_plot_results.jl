transD_GP.getchi2forall(optn)
gcf()
##
opt.xall[:] .= zall
transD_GP.plot_posterior(tempest, opt, burninfrac=0.5, figsize=(4,4),
    cmappdf="inferno")
ax = gcf().axes
ax[1].step(log10.(ρ[2:end]), z[2:end], color="w", alpha=0.4)
ax[1].invert_xaxis()
## nuisance histograms
transD_GP.plot_posterior(tempest, optn, burninfrac=0.5, nbins=50)
gcf()
