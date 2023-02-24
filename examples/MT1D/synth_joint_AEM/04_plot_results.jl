transD_GP.getchi2forall(opt)
ax = gcf().axes;
χ² = aem.ndata + length(freqs)
ax[2].plot(xlim(), [χ²/2 , χ²/2], "--", color="gray")
opt.xall[:] .= zall
transD_GP.plot_posterior(aem, opt, burninfrac=0.5, figsize=(5,6), 
    qp1=0.05, qp2=0.95, nbins=50, vmaxpc=1.0)
ax = gcf().axes
ax[1].invert_xaxis()
ax[1].step(log10.(ρ[2:end]), z[2:end], color="k", linewidth=3)
ax[1].step(log10.(ρ[2:end]), z[2:end], color="y", linewidth=1.5)
# ax[1].set_ylim(280,0)
