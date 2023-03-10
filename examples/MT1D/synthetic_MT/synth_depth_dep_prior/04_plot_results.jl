close("all")
transD_GP.getchi2forall(opt, alpha=0.5)
ax = gcf().axes;
χ² = 2*length(F.d_log10_ρ)
ax[2].plot(xlim(), [χ², χ²]/2, "--", color="gray")
## plot posterior
opt.xall[:] .= zall
transD_GP.plot_posterior(F, opt, burninfrac=0.5, qp1=0.05, qp2=0.95, nbins=50, fsize=11)
ax = gcf().axes[1]
ax.step(log10.(ρ[2:end]), z[2:end])
if F.stretch
    transD_GP.MT1DInversion.plotpriorenv(F, ax=ax, lc = "r", plotlinear=false)
end   
## plot a few models and forwards
m = transD_GP.CommonToAll.assembleTat1(opt, :fstar, temperaturenum=1)
transD_GP.MT1DInversion.plotmodelfield!(F, m[1:40:end], lcolor="k", modelalpha=0.25, logscaledepth=false)

