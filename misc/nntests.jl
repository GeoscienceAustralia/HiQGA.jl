using NearestNeighbors, Random, PyPlot, BenchmarkTools
Random.seed!(12)
# make a function to populate with nearest neighbour values
function getvoronoi!(cellvalues, kdtree, gridpoints, gridvalues)
    idxs, = knn(kdtree, gridpoints, 1)
    for i in 1:length(gridvalues)
        gridvalues[i] = cellvalues[idxs[i][1]]
    end
end
## run one 2D example
nrows, ncols = 100, 100
λ_rows, λ_cols = 1, 10.
ncells = 20
xx, yy = [i for i = 1:nrows, j = 1:ncols], [j for i = 1:nrows, j = 1:ncols]
celllocations = rand(2, ncells).*[nrows/λ_rows, ncols/λ_cols]
cellvalues = rand(ncells, 1)
kdtree = KDTree(celllocations)
gridpoints = [xx[:]/λ_rows yy[:]/λ_cols]'
gridvalues = zeros(size(xx))
getvoronoi!(cellvalues, kdtree, gridpoints, gridvalues)
figure()
imshow(gridvalues, extent=[1,ncols,nrows,1])
scatter(celllocations[2,:]*λ_cols, celllocations[1,:]*λ_rows, color="k")
xlim(extrema(xx))
ylim(extrema(yy))
## now try 3D
nrows, ncols, nnext = 80, 90, 100
λ_rows, λ_cols, λ_next = 50, 10., 2
ncells = 50
xxx, yyy, zzz = [i for i = 1:nrows, j = 1:ncols, k = 1:nnext],
                [j for i = 1:nrows, j = 1:ncols, k = 1:nnext],
                [k for i = 1:nrows, j = 1:ncols, k = 1:nnext]
celllocations = rand(3, ncells).*[nrows/λ_rows, ncols/λ_cols, nnext/λ_next]
cellvalues = rand(ncells, 1)
kdtree = KDTree(celllocations)
gridpoints = [xxx[:]/λ_rows yyy[:]/λ_cols zzz[:]/λ_next]'
gridvalues = zeros(size(xxx))
getvoronoi!(cellvalues, kdtree, gridpoints, gridvalues)
## make slices
ix, iy, iz = 50, 50, 50
fig = figure(figsize=(10,3))
ax = fig.add_subplot(1, 3, 1, projection="3d")
x, y, z = 1:nrows, 1:ncols, 1:nnext
yy, zz = [j for j = 1:ncols, k = 1:nnext],
        [k for j = 1:ncols, k = 1:nnext]
c=(gridvalues[ix,:,:].-minimum(gridvalues))./(maximum(gridvalues)-minimum(gridvalues))
ax.plot_surface(x[ix]*ones(size(yy)), yy, zz, facecolors=plt.cm.jet_r(c), shade=false)

ax = fig.add_subplot(1, 3, 2, projection="3d")
xx, zz = [i for i = 1:nrows, k = 1:nnext],
         [k for i = 1:nrows, k = 1:nnext]
c=(gridvalues[:,iy,:].-minimum(gridvalues))./(maximum(gridvalues)-minimum(gridvalues))
ax.plot_surface(xx, y[iy]*ones(size(zz)), zz, facecolors=plt.cm.jet_r(c), shade=false)

xx, yy = [i for i = 1:nrows, j = 1:ncols],
         [j for i = 1:nrows, j = 1:ncols]
ax = fig.add_subplot(1, 3, 3, projection="3d")
c=(gridvalues[:,:,iz].-minimum(gridvalues))./(maximum(gridvalues)-minimum(gridvalues))
ax.plot_surface(xx, yy, z[iz]*ones(size(yy)), facecolors=plt.cm.jet_r(c), shade=false)

for a in gcf().axes

        # a.scatter3D(celllocations[1,:][:]*λ_rows,
        #             celllocations[2,:][:]*λ_cols,
        #             celllocations[3,:][:]*λ_next,
        #             c=cellvalues[:], vmin=minimum(gridvalues), vmax=maximum(gridvalues),
        #             s=50, cmap="jet_r")
    a.set_xlim(extrema(x))
    a.set_ylim(extrema(y))
    a.set_zlim(extrema(z))
    a.set_xlabel("x")
    a.set_ylabel("y")
    a.set_zlabel("z")
    a.invert_zaxis()
end
