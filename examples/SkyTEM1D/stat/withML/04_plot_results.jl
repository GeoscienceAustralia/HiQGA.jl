close("all")
transD_GP.getchi2forall(opt)
ax = gcf().axes;
χ² = aem.ndatalow/2*log(aem.ndatalow) + aem.ndatahigh/2*log(aem.ndatahigh)
ax[2].plot(xlim(), [χ², χ²], "--", color="gray")
opt.xall[:] .= zall
transD_GP.plot_posterior(aem, opt, burninfrac=0.5, figsize=(10,6),
    cmappdf="cividis", qp1=0.01, qp2=0.99, nbins=50, vmaxpc=1.0)
ax = gcf().axes
ax[1].invert_xaxis()
ax[1].step(log10.(ρ[2:end]), z[2:end], color="r", alpha=0.5)
ax[1].step(log10.(ρ[2:end]), z[2:end], color="m", alpha=0.5, linestyle="--")
ax[2].set_ylim(280,0)
