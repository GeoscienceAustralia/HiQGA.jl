close("all")
transD_GP.getchi2forall(opt)
ax = gcf().axes;
χ² = 2*length(F.d_log10_ρ)
ax[2].plot(xlim(), [χ², χ²]/2, "--", color="gray")
##
opt.xall[:] .= zall
transD_GP.plot_posterior(F, opt, burninfrac=0.5, qp1=0.05, qp2=0.95, nbins=50, vmaxpc=1.0)
ax = gcf().axes
ax[1].step(log10.(ρ[2:end]), z[2:end], color="w", linewidth=2)
ax[1].step(log10.(ρ[2:end]), z[2:end], color="g", linestyle="--")
## plot a few models
m = transD_GP.CommonToAll.assembleTat1(opt, :fstar, temperaturenum=1)
transD_GP.MT1DInversion.plotmodelfield!(F, m[end-10:end], lcolor="k", modelalpha=0.2, logscaledepth=false)
