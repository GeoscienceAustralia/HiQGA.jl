## plot
transD_GP.getchi2forall(opt)
transD_GP.getchi2forall(optlog10λ)
transD_GP.plot_posterior(line, opt, optlog10λ, burninfrac=0.5,)
ax=gcf().axes
p = ax[1].scatter(ynoisy, x, c="w", alpha=0.2, s=25)
ax[1].plot(y, x, "--w", alpha=0.5)
del = fbounds[2]-fbounds[1]
ax[1].set_xlim(fbounds[1]-0.05del, fbounds[2]+0.05del,)
savefig("jump1D_ns.png", dpi=300)