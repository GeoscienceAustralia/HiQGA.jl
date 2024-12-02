using HiQGA
linewanted = 2008001
dr = 100
xl, yl = [], [-395,210]
vmin, vmax = -2.5, 0.5
linesmoothδ² = 1e-3
donn = false
delbin = 600
## Get LEI
line, X, Y, Z, thick, σ = transD_GP.readcols([5,7,8,9,[53,82],[23,52]], 
    "/Users/anray/Documents/work/projects/largeaem/GA_LEI/vtk_AusAEM_01_QLD/galeisbs_vector_sum_point_located_QLD.dat"; 
    startfrom=1, decfactor=5)
idx = line .== linewanted
Xd, Yd, Zd, thick, σd = map(x->x[idx,:],(X, Y, Z, thick, σ))
zalld = transD_GP.zboundarytocenter_inexact(cumsum(thick[1,:]), fudgelast=true)
## Get probabilistic
Xp, Yp, Zp, zallp, ρlow, ρmid, ρhigh, = transD_GP.readxyzrhoϕ(linewanted, 52, 
     pathname="/Users/anray/Documents/work/projects/largeaem/final_01/summaries_all/AusAEM_01/2017_QLD/summary")
## put them togather
X = [Xd, Matrix(Matrix(Xp')'), Matrix(Matrix(Xp')'), Matrix(Matrix(Xp')')]
Y = [Yd, Matrix(Matrix(Yp')'), Matrix(Matrix(Yp')'), Matrix(Matrix(Yp')')]
Z = [Zd, Matrix(Matrix(Zp')'), Matrix(Matrix(Zp')'), Matrix(Matrix(Zp')') ]
zall = [zalld, zallp, zallp, zallp]
σ = ([Matrix(log10.(σd')), -ρhigh, -ρmid, -ρlow])
## plot
titles = ["Deterministic conductivity", "10th Percentile conductivity", "50th Percentile conductivity", "90th Percentile conductivity"]
xrd, yrd, axd = transD_GP.plotmanygrids(deepcopy(σ), deepcopy(X), deepcopy(Y), deepcopy(Z), deepcopy(zall); 
    dr, vmin, vmax, donn, xl, yl, δ²=linesmoothδ², figsize=(20,10), titles, 
    preferEright=true, plotbinning=false, fontsize=10, delbin, hspace=0.22, spacefactor=0.1)
savefig("Line_$linewanted.png", dpi=500)