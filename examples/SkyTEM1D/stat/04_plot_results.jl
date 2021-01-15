close("all")
transD_GP.getchi2forall(opt)
ax = gcf().axes;
ax[2].set_ylim(10, 30)
χ² = aem.ndatalow + aem.ndatahigh
ax[2].plot(xlim(), [χ²/2 , χ²/2], "--", color="gray")
opt.xall[:] .= zall
transD_GP.plot_posterior(aem, opt, burninfrac=0.5, figsize=(4,4),
    cmappdf="inferno", qp1=0.01, qp2=0.99)
step(log10.(ρ[2:end]), z[2:end], color="w", alpha=0.4)
