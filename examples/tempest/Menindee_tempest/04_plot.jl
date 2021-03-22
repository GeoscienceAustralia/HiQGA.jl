## plot n random soundings and a background response
using  Random, PyPlot, Statistics
idx = 1
## make a closure to plot posteriors
zall, znall, zboundaries = transD_GP.CommonToAll.setupz(zstart, extendfrac, dz=dz, n=nlayers)
function plotposts(idx; computeforward=false, plotposterior=true, nbins=50,
                    computequants=false, nforwards=50, burninfrac=0.5, quants=[0.1,0.5,0.9])
    aem, znall = transD_GP.TEMPEST1DInversion.makeoperator(sounding[idx],
                               zfixed = zfixed,
                               ρfixed = ρfixed,
                               zstart = zstart,
                               extendfrac = extendfrac,
                               dz = dz,
                               ρbg = ρbg,
                               nlayers = nlayers,
                               ntimesperdecade = ntimesperdecade,
                               nfreqsperdecade = nfreqsperdecade,
                               showgeomplot = false,
                               plotfield = false)

    opt, optn = transD_GP.TEMPEST1DInversion.make_tdgp_opt(sounding[idx],
                                znall = znall,
                                fileprefix = sounding[idx].sounding_string,
                                nmin = nmin,
                                nmax = nmax,
                                K = K,
                                demean = demean,
                                sdpos = sdpos,
                                sdprop = sdprop,
                                fbounds = fbounds,
                                save_freq = save_freq,
                                λ = λ,
                                δ = δ,
                                nuisance_bounds = nuisance_bounds,
                                nuisance_sdev = nuisance_sdev,
                                updatenuisances = updatenuisances,
                                dispstatstoscreen = false
                                )
    opt.xall[:] .= zall
    if plotposterior
        transD_GP.getchi2forall(opt, alpha=0.8)
        ax = gcf().axes
        ax[2].set_ylim(10,40)
        ## plot posterior
        transD_GP.plot_posterior(aem, opt, burninfrac=burninfrac, nbins=nbins, figsize=(12,6))
        ax = gcf().axes
        ax[1].invert_xaxis()
        transD_GP.plot_posterior(aem, optn, burninfrac=burninfrac, nbins=nbins)
        ## plot forward response
    end
    if (computeforward || computequants)
        m = transD_GP.assembleTat1(opt, :fstar, temperaturenum=1, burninfrac=burninfrac)
        mn = transD_GP.CommonToAll.assemblenuisancesatT(optn, temperaturenum=1, burninfrac=burninfrac)
    end
    if computeforward
        Random.seed!(10)
        transD_GP.TEMPEST1DInversion.plotmodelfield!(aem, m[randperm(length(m))[1:nforwards]],
                                                          mn[randperm(length(m))[1:nforwards],:],
                                                          dz=dz, extendfrac=extendfrac)
    end
    if computequants
        M = vcat(m...)
        qs = hcat([quantile(M[:,i],quants) for i in 1:size(M,2)]...)
        qsn = hcat([quantile(mn[:,i],quants) for i in optn.idxnotzero]...)
        return qs, qsn
    end
end
plotposts(idx, computeforward=true)
## now try to do this for all soundings ...
quants = [0.1,0.5,0.9]
soundingquants = [plotposts(i, plotposterior=false, computequants=true, quants=quants) for i in 1:length(sounding)]
## plot them
idx = round(Int, length(sounding)/2)
X, zrx_given, xrx_given, pitchrx_given = [], [], [], []
for s in sounding
    global X, zrx_given, xrx_given, pitchrx_given
    push!(X, s.X)
    push!(zrx_given, s.z_rx)
    push!(xrx_given, s.x_rx)
    push!(pitchrx_given, s.pitch_rx)
end
X0 = X[1] - diff(X)[1]/2
xaxis = [X0, (X[1:end-1] .+ diff(X)/2)...]
cmap, vmin, vmax = "jet", -2., log10(2)
fig, ax = plt.subplots(3,1,sharex=true, sharey=true, figsize=(10,5))
a = [s[1][2,:] for s in soundingquants]
b = hcat(a...)
ax[2].pcolormesh(xaxis, zboundaries, -b, cmap=cmap, vmin=vmin, vmax=vmax)
ax[2].set_ylabel("depth m")
ax[2].set_title("P $(round((1-quants[2])*100)) conductivity")
ax[2].plot(sounding[idx].X*[1,1], ax[3].get_ylim(), "-w")
ax[2].plot(sounding[idx].X*[1,1], ax[3].get_ylim(), "--k")
a = [s[1][1,:] for s in soundingquants]
b = hcat(a...)
ax[1].pcolormesh(xaxis, zboundaries, -b, cmap=cmap, vmin=vmin, vmax=vmax)
ax[1].set_title("P $(round((1-quants[1])*100)) conductivity")
ax[1].plot(sounding[idx].X*[1,1], ax[3].get_ylim(), "-w")
ax[1].plot(sounding[idx].X*[1,1], ax[3].get_ylim(), "--k")
a = [s[1][3,:] for s in soundingquants]
b = hcat(a...)
i = ax[3].pcolormesh(xaxis, zboundaries, -b, cmap=cmap, vmin=vmin, vmax=vmax)
ax[3].plot(sounding[idx].X*[1,1], ax[3].get_ylim(), "-w")
ax[3].plot(sounding[idx].X*[1,1], ax[3].get_ylim(), "--k")
ax[3].set_title("P $(round((1-quants[3])*100)) conductivity")
xlabel("Easting m")
ax[3].set_ylim(300, 0)
transD_GP.CommonToAll.nicenup(fig, fsize=12)
fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.92, 0.15, 0.025, 0.75])
cb = fig.colorbar(i, cax=cbar_ax)
cb.ax.set_xlabel(L"\log_{10}\sigma", fontsize=12)
cb.ax.tick_params(labelsize=12)
# savefig("postice.png", dpi=300)
plotposts(idx, computeforward=true)
## median and alpha values, WIP for mesh
pl = [s[1][1,:] for s in soundingquants]
pm = [s[1][2,:] for s in soundingquants]
ph = [s[1][3,:] for s in soundingquants]
pl = hcat(pl...)
pm = hcat(pm...)
ph = hcat(ph...)
deltaP = abs(maximum(ph)-minimum(pl))
CI = abs.(ph - pl) # will be different at each pmap depth
frac = CI/deltaP
frac_max, frac_min = 0.6, 0.1 # threshholds to set
α = 0 .+ (frac_max .- frac)/(frac_max - frac_min)
idxshowall = frac .< frac_min
idxblank = frac .> frac_max
α[idxshowall] .= 1.0
α[idxblank] .= 0.
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
f = figure(figsize=(9,2))
img[zz .>topo'] .= NaN
# first fake mappable
i = imshow(img, alpha=0.5, cmap="jet", extent=[gridx[1], gridx[end], gridz[end], gridz[1]], aspect="auto")
# then blank image
imshow(ones(size(img)), extent=[gridx[1], gridx[end], gridz[end], gridz[1]], aspect="auto", cmap="gray_r")
# then actual image
imshow(img, alpha=al, cmap="jet", extent=[gridx[1], gridx[end], gridz[end], gridz[1]],aspect="auto" )
plot(X, topo, "-k")
xlabel("Easting m")
ylabel("mAHD")
ylim(-200, 152)
cb = colorbar(mappable=i)
cb.ax.set_xlabel(L"\log_{10}\sigma", fontsize=12)
cb.ax.tick_params(labelsize=12)
transD_GP.CommonToAll.nicenup(f, fsize=12)
## now try for the geometries
p1, p2, p3 = [], [], []
for s in soundingquants
    global p1, p2, p3
    push!(p1, s[2][1,:])
    push!(p2, s[2][2,:])
    push!(p3, s[2][3,:])
end
p1 = hcat(p1...)
p2 = hcat(p2...)
p3 = hcat(p3...)
f, ax = plt.subplots(3,1, sharex=true, figsize=(10,2.6))
ax[1].plot(X, p1[1,:])
ax[1].plot(X, p2[1,:])
ax[1].plot(X, p3[1,:])
ax[1].plot(X, zrx_given, "--k")
ax[1].set_ylabel("height m")
ax[1].set_title("zrx")
ax[2].plot(X, p1[2,:])
ax[2].plot(X, p2[2,:])
ax[2].plot(X, p3[2,:])
ax[2].plot(X, xrx_given, "--k")
ax[2].set_ylabel("tx-rx sep m")
ax[2].set_title("xrx")
ax[3].plot(X, p1[3,:])
ax[3].plot(X, p2[3,:])
ax[3].plot(X, p3[3,:])
ax[3].plot(X, pitchrx_given, "--k")
ax[3].set_ylabel("degrees")
ax[3].set_title("pitch_rx")
ax[3].set_xlabel("Easting m")
# ax[3].invert_xaxis()
transD_GP.nicenup(f)
