using NearestNeighbors
xx, zz = [x for z in zall, x in X], [z for z in zall, x in X]
topo = [s.Z for s in sounding]
zz = topo' .- zz # mAHD
kdtree = KDTree([xx[:]'; zz[:]'])
gridx = range(X[1], X[end], length=length(X))
gridz = reverse(range(extrema(zz)..., step=dz))
xx, zz = [x for z in gridz, x in gridx], [z for z in gridz, x in gridx]
idxs, = nn(kdtree, [xx[:]'; zz[:]'])
img = zeros(size(xx))
al = zeros(size(xx))
for i = 1:length(img)
    img[i] = -pm[idxs[i]]
    al[i] = α[idxs[i]]
end
f = figure(figsize=(8,5))
img[zz .>topo'] .= NaN
imshow(img, alpha=al, cmap="jet", extent=[gridx[1], gridx[end], gridz[end], gridz[1]], aspect="auto")
plot(X, topo, "-k")
