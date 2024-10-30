dr = 6.25
lnames = [100401, 100502]
idx = [[200, 900], [60]]
burninfrac = 0.25
qp1, qp2 = 0.1, 0.9
vmin, vmax = -reverse(fbounds[:])
transD_GP.summaryAEMimages(soundings, opt;yl=[-150,100],
    zall, dr, dz, cmap="turbo", burninfrac, figsize=(12,8), fontsize=18,
    vmin=vmin, vmax=vmax, lnames, idx, useML, qp1, qp2, dpi=400,
    preferEright=true, saveplot=false)
## a few individual soundings
transD_GP.plotindividualAEMsoundings(soundings, aem, opt, idx;
    burninfrac, nbins=50, computeforwards=true, zall, qp1, qp2, 
    nforwards=10, usekde=true, lnames)