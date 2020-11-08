transD_GP.getchi2forall(opt, nxticks=4, fsize=6)
ax = gcf().axes;
ax[2].set_ylim(10, 30)
transD_GP.getchi2forall(optlog10λ, nxticks=4, fsize=6)
ax = gcf().axes;
ax[2].set_ylim(10, 30)
opt.xall[:] .= zall
transD_GP.plot_posterior(aem, opt, optlog10λ, burninfrac=0.5,
    figsize=(7.8,4), cmappdf="inferno", qp1=0.01, qp2=0.99)
ax = gcf().axes
ax[1].step(log10.(ρ[2:end]), z[2:end], alpha=0.4, color="w")
