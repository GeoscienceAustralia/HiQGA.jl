close("all")
transD_GP.getchi2forall(opt)
ax = gcf().axes;
ax[2].set_ylim(10, 30)
χ² = aem.ndatalow + aem.ndatahigh
ax[2].plot(xlim(), [χ²/2 , χ²/2], "--", color="gray")
opt.xall[:] .= zall
transD_GP.plot_posterior(aem, opt, burninfrac=0.5, figsize=(10,6),
    cmappdf="inferno", qp1=0.05, qp2=0.95, nbins=200, vmaxpc=0.5)
ax = gcf().axes
ax[1].invert_xaxis()
ax[1].step(log10.(ρ[2:end]), z[2:end], color="w", alpha=0.8)
