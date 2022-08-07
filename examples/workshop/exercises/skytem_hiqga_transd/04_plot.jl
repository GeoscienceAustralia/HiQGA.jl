dxplot     = 100
idx        = [10]
burninfrac = 0.25
vmin, vmax = -reverse([extrema(fbounds)...])
transD_GP.SkyTEM1DInversion.summaryimages(soundings, opt,
    zstart        = zstart,
    d             = dz,
    dr            = dxplot,
    cmap          = "jet",
    extendfrac    = extendfrac,
    vmin          = vmin, 
    vmax          = vmax,
    idx           = idx,
    nlayers       = nlayers,
    preferEright  = true,
    burninfrac    = burninfrac,
    saveplot      = true,
    figsize       = (12,12),
    yl            = [-50,100]
    )
##
computeforwards = true
nforwards       = 1
nbins           = 50
transD_GP.SkyTEM1DInversion.plotindividualsoundings(soundings, opt,
    burninfrac=burninfrac,
    nbins           = nbins,
    zfixed          = zfixed,
    ρfixed          = ρfixed,
    zstart          = zstart,
    extendfrac      = extendfrac,
    dz              = dz,
    ρbg             = ρbg,
    nlayers         = nlayers,
    computeforwards = computeforwards,
    nforwards       = nforwards,
    idxcompute      = idx,
    pdfclim         = [0, 2.5])