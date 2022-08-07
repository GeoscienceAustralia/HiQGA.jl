dr = 15
idx = [415]
burninfrac=0.25
vmin, vmax = -2.5, 0.5
transD_GP.SkyTEM1DInversion.summaryimages(soundings, opt,
                   zstart=zstart,
                   dz=dz,
                   dr=dr,
                   cmap="jet",
                   extendfrac=extendfrac,
                   vmin=vmin, vmax=vmax,
                   idx=idx,
                   nlayers=nlayers,
                   useML=useML,
                   preferEright=true,
                   burninfrac = burninfrac,
                   saveplot=true,
                   topowidth = 1,
                   figsize=(12,12),
                   omitconvergence = false,
                   yl = [-50,100]
                   )