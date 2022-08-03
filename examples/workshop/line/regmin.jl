using LinearAlgebra, PyPlot, HiQGA, LinearMaps, Random
## let's use this for a linear problem to prove a point
Random.seed!(11)
x = -5:.1:5
nuse = 7
yclean = 0.5(sin.(x).^2 .+ cos.(x).^2)
useidx = randperm(length(yclean))[1:nuse]
d = (yclean + 0.1*randn(length(x)))[useidx] 
Gfull = [sin.(x).^2 cos.(x).^2]
G = Gfull[useidx,:]
R = LinearMap(transD_GP.R1Dop, transD_GP.Rt1Dop, 2)
m = [Matrix(G'G + δ²*R'R)\(G'd) for δ² in 10 .^[-16., 0, 12, 16]]
d_pred = [Gfull*m_ for m_ in m]
##
figure()
subplot(121)
plot(x, yclean, "--k", label="true")
plot(x[useidx][sortperm(x[useidx])], d[sortperm(x[useidx])], ".-", label = "noisy")
[plot(x, d_pred[i], label="$(m[i])") for i in 1:length(m)]
# legend()
# for error ellipses
x1 = 0.3:.01:0.6
x2 = 0.3:.01:0.6
s2 = subplot(122)
plot(x1, x2, "-k")
mall = reduce(hcat, m)'
plot(mall[:,1], mall[:,2], ".-")
s2.set_aspect(1.)
s2.axis([extrema(x1)...,extrema(x2)...])
xx2, xx1 = [x for x in x1, y in x2],[y for x in x1, y in x2] 
rmat = xx1.^2 + xx2.^2 - 2*xx1.*xx2
##
figure(figsize=(16,6))
d_pred = G*[xx1[:]';xx2[:]']
res = reshape(sum((d_pred .- d).^2, dims=1), size(xx1))
s1 = subplot(131, sharex=s1, sharey=s1)
imshow(res, extent=[x1[1],x1[end],x2[1],x2[end]], aspect="equal", origin="lower", cmap="jet")
s1.plot(x1,x2, "-k", linewidth=1)
s1.plot(mall[1,1], mall[1,2], "+w", markersize=15)
s1.plot(mall[1,1], mall[1,2], ".w", markersize=15)
xlabel("x1"); ylabel("x2")
s2 = subplot(132, sharex=s1, sharey=s1)
imshow(rmat, extent=[x1[1],x1[end],x2[1],x2[end]], aspect="equal", origin="lower", cmap="hot")
s2.plot(x1,x2, "-w", linewidth=1)
s2.plot(mall[1,1], mall[1,2], "+w", markersize=15)
s2.plot(mall[1,1], mall[1,2], ".w", markersize=15)
xlabel("x1"); ylabel("x2")
s3 = subplot(133, sharex=s1, sharey=s1)
imshow(res, extent=[x1[1],x1[end],x2[1],x2[end]], aspect="equal", origin="lower", cmap="jet", alpha=0)
xlabel("x1"); ylabel("x2")
cp = contour(xx1, xx2, rmat, linewidth=2.0, cmap="hot")
cp = contour(xx1, xx2, res, linewidth=2.0, cmap="jet")
s3.plot(x1,x2, "-k", linewidth=1)
plot(mall[1:end-1,1], mall[1:end-1,2], ".--k", markersize=15)
s3.plot(mall[1,1], mall[1,2], "+w", markersize=15)
s1.set_xticklabels([]); s1.set_yticklabels([])
transD_GP.nicenup(gcf(), fsize=24)
# s3.clabel(cp, inline=1, fontsize=10)
