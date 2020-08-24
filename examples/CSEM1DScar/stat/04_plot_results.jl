GeophysOperator.getchi2forall(opt, fsize=8)
ax = gcf().axes;
ax[3].set_xlim(2000001, 4000001)
ax[2].set_ylim(2050, 2310)
plt.tight_layout()
savefig("csem_scar_conv_s.png", dpi=300)
opt.xall[:] .= zall
GeophysOperator.plot_posterior(csem, opt, burninfrac=0.5, fsize=8,
    figsize=(4,4), cmappdf="inferno", nbins=100, qp1=0.01, qp2=0.99)
ax = gcf().axes
ax[1].step(log10.(ρ[2:end]), z[2:end], alpha=0.4, color="w")
ax[1].axis([-0.8591188164842742, 2.901547748868629, z[end-1], z[nfixed+1]])
gcf().text(0.02, 0.9, "a.", fontsize=14, color="red")
savefig("csem_scar_post_s.png", dpi=300)
M = GeophysOperator.assembleTat1(opt, :fstar, temperaturenum=1)
Random.seed!(10)
GeophysOperator.plotmodelfield!(csem, M[randperm(length(M))[1:30]],
                                       dz=dz, extendfrac=extendfrac, onesigma=false)
gcf().text(0.02, 0.9, "a.", fontsize=14, color="red")
savefig("csem_scar_fwds_post_s.png", dpi=300)
