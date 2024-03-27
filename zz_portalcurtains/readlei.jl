using HiQGA
linewanted = 2008001
dr = 100
xl, yl = [], []
vmin, vmax = -2.5, 0.5
linesmoothδ² = 1e-2
donn = false
delbin = 700
## Get LEI
line, X, Y, Z, thick, σ = transD_GP.readcols([5,7,8,9,[53,82],[23,52]], 
    "/Users/anray/Documents/work/projects/largeaem/GA_LEI/vtk_AusAEM_01_QLD/galeisbs_vector_sum_point_located_QLD.dat"; 
    startfrom=1, decfactor=48)
idx = line .== linewanted
Xd, Yd, Zd, thick, σd = map(x->x[idx,:],(X, Y, Z, thick, σ))
zalld = transD_GP.zboundarytocenter(cumsum(thick[1,:]))
## another adjacent line
linewanted = 2006001
line, X, Y, Z, thick, σ = transD_GP.readcols([5,7,8,9,[53,82],[23,52]], 
    "/Users/anray/Documents/work/projects/largeaem/GA_LEI/vtk_AusAEM_01_QLD/galeisbs_vector_sum_point_located_QLD.dat"; 
    startfrom=1, decfactor=48)
idx = line .== linewanted
Xp, Yp, Zp, thick, σp = map(x->x[idx,:],(X, Y, Z, thick, σ))
zallp = transD_GP.zboundarytocenter(cumsum(thick[1,:]))
## put them togather
X = [Xd, Matrix(Matrix(Xp')')]
Y = [Yd, Matrix(Matrix(Yp')')]
Z = [Zd, Matrix(Matrix(Zp')')]
zall = [zalld, zallp]
σ = ([Matrix(log10.(σd')), Matrix(log10.(σp'))])
## plot
xrd, yrd, axd = transD_GP.plotmanygrids(σ, X, Y, Z, zall; dr, vmin, vmax, donn, xl, yl, δ²=linesmoothδ²,
    preferEright=true, plotbinning=true, fontsize=10, delbin)