close("all")
transD_GP.getchi2forall(opt)
ax = gcf().axes;
χ² = length(F.d_log10_ρ)
ax[2].plot(xlim(), [χ², χ²], "--", color="gray")
ax[2].set_ylim(χ²-10, χ²+10)
##
opt.xall[:] .= zall
transD_GP.plot_posterior(F, opt, burninfrac=0.5, figsize=(10,6), qp1=0.05, qp2=0.95, nbins=50, vmaxpc=1.0)
ax = gcf().axes
ax[1].step(log10.(ρ[2:end]), z[2:end], color="w", linewidth=2)
ax[1].step(log10.(ρ[2:end]), z[2:end], color="g", linestyle="--")

