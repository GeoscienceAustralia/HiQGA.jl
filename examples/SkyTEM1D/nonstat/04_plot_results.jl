# convergence statistics
close("all")
transD_GP.getchi2forall(opt)
ax = gcf().axes;
ax[2].set_ylim(10, 30)
χ² = aem.ndatalow + aem.ndatahigh
ax[2].plot(xlim(), [χ²/2 , χ²/2], "--", color="gray")
transD_GP.getchi2forall(optlog10λ)
ax = gcf().axes;
ax[2].set_ylim(10, 30)
ax[2].plot(xlim(), [χ²/2 , χ²/2], "--", color="gray")
# posterior
opt.xall[:] .= zall
transD_GP.plot_posterior(aem, opt, optlog10λ, burninfrac=0.5,
    figsize=(15,6), cmappdf="inferno", nbins=100, qp1=0.05, qp2=0.95, vmaxpc=0.5)
ax = gcf().axes
ax[1].invert_xaxis()
ax[1].step(log10.(ρ[2:end]), z[2:end], color="w")
ax[1].step(log10.(ρ[2:end]), z[2:end], color="k", linestyle="--")
