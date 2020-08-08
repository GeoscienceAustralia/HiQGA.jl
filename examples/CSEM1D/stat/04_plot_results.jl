GeophysOperator.getchi2forall(opt, nxticks=4, fsize=6)
ax = gcf().axes;
ax[3].set_ylim(100, 200)
ax[4].set_ylim(100, 200)
savefig("csem_conv_s.png", dpi=300)
opt.xall[:] .= zall
GeophysOperator.plot_posterior(csem, opt, burninfrac=0.5, figsize=(4,4),
    cmappdf="inferno")
step(log10.(ρ[2:end]), z[2:end], color="w", alpha=0.6)
ax = gcf().axes
ax[1].set_xlim(-1.0959789629712384, 2.5166552528237873)
savefig("csem_post_s.png", dpi=300)
