## stats
GeophysOperator.getchi2forall(opt, alpha=0.7)
ax = gcf().axes
ax[2].set_ylim(10, 50)
## posterior models
opt.xall[:] .= zall
GeophysOperator.plot_posterior(aem, opt, burninfrac=0.5, fsize=10,
           figsize=(4,6), cmappdf="inferno", nbins=100, qp1=0.1, qp2=0.9)
M = GeophysOperator.assembleTat1(opt, :fstar, temperaturenum=1)
Random.seed!(10)
GeophysOperator.SkyTEM1DInversion.plotmodelfield!(aem, M[randperm(length(M))[1:50]],
                                       dz=dz, extendfrac=extendfrac, onesigma=false, alpha=0.05)
