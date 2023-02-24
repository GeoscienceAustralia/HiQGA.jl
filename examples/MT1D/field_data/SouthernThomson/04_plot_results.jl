using PyPlot
transD_GP.getchi2forall(opt, alpha=0.5)
ax = gcf().axes;
χ² = length(F.d_log10_ρ)
F.useML || ax[2].plot(xlim(), [χ², χ²], "--", color="gray")
F.useML || ax[2].set_ylim(χ²-10, χ²+30)
## plot posterior
opt.xall[:] .= zall
transD_GP.plot_posterior(F, opt, burninfrac=0.2, figsize=(5,6), qp1=0.1, qp2=0.9, nbins=50, fsize=11)
## plot a few models
m = transD_GP.CommonToAll.assembleTat1(opt, :fstar, temperaturenum=1)
transD_GP.MT1DInversion.plot_posterior(F, m[1:40:end], lcolor="k", modelalpha=0.2)

