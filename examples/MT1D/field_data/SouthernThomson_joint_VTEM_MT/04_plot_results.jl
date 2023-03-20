using PyPlot
## plot misfit
transD_GP.getchi2forall(opt, omittemp=true, figsize=(5,3))
savefig("aa_"*fname*"_misfit.png", dpi=300)
## plot posterior
opt.xall[:] .= zall
transD_GP.plot_posterior(aem, opt, burninfrac=0.2, figsize=(5,6), 
    qp1=0.1, qp2=0.9, nbins=50, vmaxpc=1.0, plotmean=false)
ax = gcf().axes
ax[1].invert_xaxis()
ax[1].set_title(fname)
savefig("bb_"*fname*"_posterior.png", dpi=300)
ax[1].set_ylim(286.5, 0.75)
savefig("bb_"*fname*"_posterior_shallow.png", dpi=300)
## get models and plot fit
m = transD_GP.CommonToAll.assembleTat1(opt, :fstar, temperaturenum=1)
transD_GP.VTEM1DInversion.plotmodelfield!(aem, m[end-50:end], figsize=(6,3))
savefig("cc_"*fname*"_AEM_forwards.png", dpi=300)
transD_GP.MT1DInversion.plotmodelfield!(FMT, m[end-50:end], logscaledepth=false, revax=true, figsize=(9,3))
savefig("dd_"*fname*"_MT_forwards.png", dpi=300)
## bash script put into one png file
open("pngscript.sh", "w") do file
    write(file, "#!/bin/bash\nmontage aa*png bb*png cc*png dd*png -tile 3x2 -geometry +0+0 "*fname*"_combined.png\n")
    write(file, "#!/bin/bash\nmontage aa*png bb*shallow*png cc*png dd*png -tile 2x2 -geometry +0+0 "*fname*"_combined_shallow.png\n")
end
run(`chmod +x pngscript.sh`)
run(`./pngscript.sh`)

