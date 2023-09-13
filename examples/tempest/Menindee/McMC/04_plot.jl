
idxplot = [1,2]
burninfrac = 0.5
nbins = 50
computeforwards, nforwards = true, 2
interpolatedr = 10
vmax, vmin = .- extrema(fbounds)
## get the summary images of resistivity and nuisances over the profile
transD_GP.AEMwithNuisanceMcMCInversionTools.summaryAEMwithnuisanceimages(soundings, opt, optn; zall,
                            burninfrac, dz, useML, dr = interpolatedr, fontsize = 10,
                            vmin, vmax, cmap="turbo", figsize=(8,10), idx = idxplot,
                            omitconvergence = false, preferEright = false,  preferNright = false,
                            saveplot = true,  labelnu = ["zRx m", "xRx m"], vectorsum 
                        )

## now plot individual soundings in in idxplot
transD_GP.AEMwithNuisanceMcMCInversionTools.plotindividualsoundings(soundings, aem, opt, optn, idxplot;
                        burninfrac,nbins, figsize  = (12,6), zall,
                        computeforwards, nforwards,
                        )