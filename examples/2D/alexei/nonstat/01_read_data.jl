## get data x, y, property
using DelimitedFiles, PyPlot
d, h = readdlm("mohodepth_GDA_LCC.csv",',', header=true)
x, y = d[:,5], d[:,6]
confidence = d[:,3]
confidence[confidence.<0.1].=0.1
σ = 1 ./confidence
z = d[:,4]
f = figure()
scatter(x, y, c=z, cmap="RdBu")
xlabel("Easting m")
ylabel("Northing m")
cb = colorbar()
cb.ax.set_xlabel("km")
gca().set_aspect(1.)
## make a grid
using NearestNeighbors
xmin, xmax = extrema(x)
ymin, ymax = extrema(y)
delx, dely = 20_000., 20_000
xgrid, ygrid = [x for y in ymin:dely:ymax, x in xmin:delx:xmax], [y for y in ymin:dely:ymax, x in xmin:delx:xmax]
xall = [xgrid[:]'; ygrid[:]']
kdtree = KDTree(xall)
idxs, = nn(kdtree, [x[:]';y[:]'])
imagegrid = NaN .+ zeros(size(xgrid))
noisegrid = copy(imagegrid)
imagegrid[idxs] .= z
noisegrid[idxs] .= σ
figure(figsize=(12,4))
s1 = subplot(121)
imshow(imagegrid, extent=[xmin, xmax, ymin, ymax], origin="lower", aspect=1., cmap="RdBu")
title("moho depth")
xlabel("Easting m")
ylabel("Northing m")
cb = colorbar()
cb.ax.set_xlabel("km")
subplot(122, sharex=s1, sharey=s1)
imshow(noisegrid, extent=[xmin, xmax, ymin, ymax], origin="lower", aspect=1., cmap="viridis")
title("noise on depth")
xlabel("Easting m")
ylabel("Northing m")
cb = colorbar()
cb.ax.set_xlabel("km")
plt.tight_layout()