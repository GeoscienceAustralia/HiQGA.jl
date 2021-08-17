## plot n random soundings and a background response
using  Random, PyPlot, Statistics, NearestNeighbors, DelimitedFiles
idx = 1
## make a closure to plot summary posteriors
zall, znall, zboundaries = transD_GP.CommonToAll.setupz(zstart, extendfrac, dz=dz, n=nlayers)
function summarypost(;qp1=0.05, qp2=0.95, nbins=50, burninfrac=0.5)
    aem, znall = transD_GP.TEMPEST1DInversion.makeoperator(sounding[1],
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

    opt, optn = transD_GP.TEMPEST1DInversion.make_tdgp_opt(sounding[1],
                                znall = znall,
                                fileprefix = sounding[1].sounding_string,
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
    pl, pm, ph, ρmean, vdmean, vddev = map(x->zeros(length(zall), length(sounding)), 1:6)
    for idx = 1:length(sounding)
        @info "$idx out of $(length(sounding))"
        opt.fdataname = sounding[idx].sounding_string*"_"
        pl[:,idx], pm[:,idx], ph[:,idx], ρmean[:,idx],
        vdmean[:,idx], vddev[:,idx] = transD_GP.plot_posterior(aem, opt, burninfrac=burninfrac,
                                                qp1=qp1, qp2=qp2,
                                                nbins=nbins, doplot=false)
    end
    pl, pm, ph, ρmean, vdmean, vddev
end
# interpolate grid to mAHD with topo
function makegrid(a;horizfillfactor=5)
    X = [s.X for s in sounding]
    xx, zz = [x for z in zall, x in X], [z for z in zall, x in X]
    topo = [s.Z for s in sounding]
    zz = topo' .- zz # mAHD
    kdtree = KDTree([xx[:]'; zz[:]'])
    gridx = range(X[1], X[end], length=round(Int, horizfillfactor)*length(X))
    gridz = reverse(range(extrema(zz)..., step=dz))
    xx, zz = [x for z in gridz, x in gridx], [z for z in gridz, x in gridx]
    idxs, = nn(kdtree, [xx[:]'; zz[:]'])
    img = zeros(size(xx))
    for i = 1:length(img)
        img[i] = a[idxs[i]]
    end
    kdtree = KDTree(X[:]')
    idxs, = nn(kdtree, gridx[:]')
    topofine = [topo[idxs[i]] for i = 1:length(gridx)]
    img[zz .>topofine'] .= NaN
    img, gridx, gridz, topofine
end
# this function does too many things
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
function makeα(pl, ph; frac_max=0.6, frac_min=0.1)
    deltaP = abs(maximum(ph)-minimum(pl))
    CI = abs.(ph - pl) # will be different at each pmap depth
    frac = CI/deltaP
    α = 0 .+ (frac_max .- frac)/(frac_max - frac_min)
    idxshowall = frac .< frac_min
    idxblank = frac .> frac_max
    α[idxshowall] .= 1.0
    α[idxblank] .= 0.
    α
end
function plotprofile(ax, idx)
    ax.plot(sounding[idx].X*[1,1], [ax.get_ylim()[1], sounding[idx].Z], "-w")
    ax.plot(sounding[idx].X*[1,1], [ax.get_ylim()[1], sounding[idx].Z], "--k")
end
function plotgrids(grids...; gridx=[1, 2.], gridz=[1,2.], vmin=-2, vmax=log10(2),
            elev=nothing, ylabels=nothing,
            figsize=(18,6), α=1., cmap="jet", titles=nothing, idx=nothing, fsize=12)
    f, ax = plt.subplots(size(grids, 1), 1, figsize=figsize,
                        sharex=true, sharey=true, squeeze=false)
    for (i, g) in enumerate(grids)
        # first fake mappable
        img = ax[i].imshow(g, cmap=cmap, aspect="auto", vmax=vmax, vmin = vmin,
                extent=[gridx[1], gridx[end], gridz[end], gridz[1]], alpha=0.5)
        # then blank image
        ax[i].imshow(ones(size(g)), extent=[gridx[1], gridx[end], gridz[end], gridz[1]], aspect="auto", cmap="gray_r")
        # then actual image
        ax[i].imshow(g, cmap=cmap, aspect="auto", alpha=α, vmax=vmax, vmin = vmin,
                extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
        colorbar(img, ax=ax[i])
        titles == nothing || ax[i].set_title(titles[i])
        idx  == nothing || plotprofile.(ax[i], idx)
        elev == nothing || plot(gridx, elev, "-k")
        ylabels == nothing || ax[i].set_ylabel(ylabels[i])
    end
    ax[end].set_xlabel("Easting m")
    transD_GP.CommonToAll.nicenup(f, fsize=fsize)
end
## use these functions
qp1, qp2=0.05, 0.95
pl, pm, ph, ρmean, vdmean, vddev = summarypost(burninfrac=0.5, qp1=qp1, qp2=1p2)
writedlm("plow.txt", pl)
writedlm("pmid.txt", pm)
writedlm("phigh.txt", ph)
writedlm("mean.txt", ρmean)
writedlm("vdmean.txt", vdmean)
writedlm("vddev.txt", vddev)

α = makeα(pl, ph, fracmax=0.6, fracmin=0.1)
αgrid, gridx, gridz, topofine = makegrid(α)

phgrid, = makegrid(-pl)
pmgrid, = makegrid(-pm)
plgrid, = makegrid(-ph)

α_alt = vdmean./vddev
α_alt /= maximum(α_alt)
α_alt .^=1/2.2

α_altgrid, = makegrid(α_alt)

idxs=[140, 300]
plotgrids(plgrid, pmgrid, phgrid; gridx=gridx, gridz=gridz, elev=nothing, ylabels=["mAHD", "mAHD", "mAHD"],
            α=1, titles=["Percentile $(100*qp1) conductivity", "Percentile 50 conductivity", "Percentile $(100*qp2) conductivity"],
            cmap="viridis", idx=idxs, vmin=-1.5, vmax=log10(2))
ylim(-200,100)
plotposts.(idxs)
