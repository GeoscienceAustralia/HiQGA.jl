cd(@__DIR__)
using HiQGA, PyPlot
linetoget = 1009001
include("0000_00_get_linelist.jl")
idxgood = [s.linenum == linetoget for s in soundings]
soundings = soundings[idxgood]
include("0000_02_set_options.jl")
##
lnames = [1009001]
idx = [90]
burninfrac = 0.25
nbins = 50
computeforwards, nforwards = true, 10 
interpolatedr = 100 
vmax, vmin = .- extrema(fbounds)
qp1, qp2 = 0.1, 0.9 
## get the summary images of resistivity and nuisances over the profile
transD_GP.summaryAEMimages(soundings, opt, optn; zall, qp1, qp2, lnames, idx,
                            burninfrac, dz, useML, dr = interpolatedr, fontsize = 8, yl = [], 
                            vmin, vmax, cmap="turbo", figsize=(20,5), dpi=300,
                            omitconvergence = false, preferEright = true,  preferNright = false, showplot = true,
                            saveplot = true,  labelnu = ["zRx m", "Tx-Rx hsep m"], vectorsum, Rmax=nothing, 
                        )
##
## now plot individual soundings in in idxplot
transD_GP.AEMwithNuisanceMcMCInversionTools.plotindividualsoundings(soundings, aem, opt, optn, idx;
                        burninfrac,nbins, figsize  = (6,6), zall, linecolor="k",qp1, qp2,
                        computeforwards, nforwards, alpha=0.5,
                        )
## now get something from the ASEG-GDF and plot it
transD_GP.dfn2hdr("/scratch/ns59/ar0754/AusAEM_03_ERC/AusAEM_03_ERC_2021_ERC_01_EPSG_28354.dfn")
# from DFN file
# 1    Line
# 2    X   
# 3    Y   
# 4    Z   
# 5   56   zcenter
# 57  108  log10_cond_low
# 109 160  log10_cond_mid
# 161 212  log10_cond_high
# 213 264  log10_cond_avg
# 265  phid_mean
# 266  phid_sdev
# 267  z_rx_low
# 268  x_rx_low
# 269  z_rx_mid
# 270  x_rx_mid
# 271  z_rx_high
# 272  x_rx_high
A = transD_GP.readlargetextmatrix("/scratch/ns59/ar0754/AusAEM_03_ERC/AusAEM_03_ERC_2021_ERC_01_EPSG_28354.dat")
idxgood = A[:,1] .== linetoget
A = A[idxgood,:]
mhigh = -A[idx, 57:108][:] # res not cond
mmid  = -A[idx,109:160][:]
mlow  = -A[idx,161:212][:]
aemuse = makeoperator(aem, soundings[idx][1])
nuinv = [A[idx,267:268][:], A[idx,269:270][:], A[idx,271:272][:]]
# nugiven = transD_GP.TEMPEST1DInversion.getnufromsounding(soundings[idx[1]])
mn = reduce(hcat, [transD_GP.TEMPEST1DInversion.setnuforinvtype(aemuse, n) for n in nuinv])'
transD_GP.plotmodelfield!(aemuse, [mlow,mmid,mhigh], collect(mn); model_lw=1.5, forward_lw=1.5, alpha=1)
# misfit for median model
transD_GP.TEMPEST1DInversion.getfield!(mmid, mn[2,:], aemuse)
phid = transD_GP.TEMPEST1DInversion.getchi2by2(aemuse)/(aem.ndatax/2)
@info "phid is $phid"
