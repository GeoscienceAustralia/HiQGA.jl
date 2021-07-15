close("all")
transD_GP.getchi2forall(opt, fsize=6)
ax = gcf().axes;
ax[2].set_ylim(100, 200)
χ² = sum(2*csem.ndata)
ax[2].plot(xlim(), [χ²/2 , χ²/2], "--", color="gray")
savefig("csem_conv_s.png", dpi=300)
opt.xall[:] .= zall
transD_GP.plot_posterior(csem, opt, burninfrac=0.5, figsize=(4,4),
    cmappdf="inferno", qp1=0.1, qp2=0.9, nbins=100)
ax = gcf().axes
ax[1].step(log10.(ρ[2:end]), z[2:end], color="w", alpha=0.4)
ax[1].set_xlim(-1.0959789629712384, 2.5166552528237873)
savefig("csem_post_s.png", dpi=300)
