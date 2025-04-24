cd(@__DIR__)
# ENV["MPLBACKEND"]="agg"
using HiQGA, PyPlot, Printf
linewanted = 4019001 #4014001, figsize=(25,6), yl=[-280, 170]
dr = 100
xl, yl = [], [-180, 220]
vmin, vmax = -2.5, 0.5
linesmoothδ² = 1e-3
donn = false
delbin = 500
## Get LEI
line, X, Y, Z, thick, σ = transD_GP.readcols([5,7,8,9,[53,82],[23,52]], 
    "/scratch/ns59/ar0754/keep/All_AEM/AusAEM20_Tempest/Earaheedy& DesertStrip_AusAEM2020_GA_vsum_inversion.dat",
    startfrom=1, decfactor=5)
idx = line .== linewanted
Xd, Yd, Zd, thick, σd = map(x->x[idx,:],(X, Y, Z, thick, σ))
zalld = transD_GP.thicktodepth(thick[1,:])
## Get probabilistic
Xp, Yp, Zp, zallp, ρlow, ρmid, ρhigh, = transD_GP.readxyzrhoϕ(linewanted, 52, 
     pathname="/g/data/z67/ar0754/largeaem/production/AusAEM20_WA/Earaheedy_Desert_Strip/summary/")
## put them together
X = [Xd, Matrix(Matrix(Xp')'), Matrix(Matrix(Xp')'), Matrix(Matrix(Xp')')]
Y = [Yd, Matrix(Matrix(Yp')'), Matrix(Matrix(Yp')'), Matrix(Matrix(Yp')')]
Z = [Zd, Matrix(Matrix(Zp')'), Matrix(Matrix(Zp')'), Matrix(Matrix(Zp')') ]
zall = [zalld, zallp, zallp, zallp]
σ = ([Matrix(log10.(σd')), -ρhigh, -ρmid, -ρlow])
## plot
XYprofiles = nothing#[Xp[520];Yp[520]]
titles = ["Deterministic conductivity", "10th Percentile conductivity", "50th Percentile conductivity", "90th Percentile conductivity"]
xrd, yrd, axd = transD_GP.plotmanygrids(deepcopy(σ), deepcopy(X), deepcopy(Y), deepcopy(Z), deepcopy(zall); 
    dr, vmin, vmax, donn, xl, yl, δ²=linesmoothδ², figsize=(12.5,6), titles, 
    preferEright=true, plotbinning=false, fontsize=11, delbin, hspace=0.6, smallratio=.3, spacefactor=.015)
VE = @sprintf("%.1f", transD_GP.getVE(axd[1]))
savefig("Line_$(linewanted)_VE_$VE.png", dpi=500)
