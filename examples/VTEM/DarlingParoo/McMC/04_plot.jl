dr = 20
idx = [1,2]
burninfrac = 0.5
vmin, vmax = -2.5, 0.5
transD_GP.summaryAEMimages(soundings, opt;
    zall, dr, dz, cmap="turbo", burninfrac,
    vmin=vmin, vmax=vmax, idx, useML,
    preferEright=true, saveplot=true)
##
computeforwards = true
nforwards = 5
nbins = 50
transD_GP.plotindividualAEMsoundings(soundings, aem, opt, idx;
    burninfrac, nbins, computeforwards, zall,
    nforwards)

