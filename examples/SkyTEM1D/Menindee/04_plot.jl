dr = 20
idx = 1:100:length(soundings)
burninfrac=0.5
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
                   saveplot=true
                   )

computeforwards = true
nforwards = 100
nbins = 50
transD_GP.SkyTEM1DInversion.plotindividualsoundings(soundings, opt,
    burninfrac=burninfrac,
    nbins = nbins,
    zfixed = zfixed,
    ρfixed = ρfixed,
    zstart = zstart,
    extendfrac = extendfrac,
    dz = dz,
    ρbg = ρbg,
    nlayers = nlayers,
    ntimesperdecade = ntimesperdecade,
    nfreqsperdecade = nfreqsperdecade,
    computeforwards = computeforwards,
    nforwards = nforwards,
    idxcompute = idx)

