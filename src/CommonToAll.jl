module CommonToAll
using TransD_GP, PyPlot, StatsBase, Statistics,
      AbstractOperator

export trimxft, assembleTat1, gettargtemps, checkns, getchi2forall, nicenup,
        plot_posterior, make1Dhist, make1Dhists, setupz, plotdepthtransforms,
        unwrap, getn, geomprogdepth

function trimxft(opt::TransD_GP.Options, burninfrac::Float64)
    x_ft = assembleTat1(opt, :x_ftrain, burninfrac=burninfrac)
    n = assembleTat1(opt, :nodes, burninfrac=burninfrac)
    x, ft = zeros(size(opt.xall, 1), sum(n)), zeros(size(opt.fbounds, 1), sum(n))
    nlast = 0
    for (i, xft) in enumerate(x_ft)
        x[:,nlast+1:nlast+n[i]]  = xft[1:size(opt.xall, 1), 1:n[i]]
        ft[:,nlast+1:nlast+n[i]] = xft[size(opt.xall, 1)+1:end, 1:n[i]]
        nlast += n[i]
    end
    x, ft
end

function assembleTat1(optin::TransD_GP.Options, stat::Symbol; burninfrac=0.5)
    isns = checkns(optin)
    @assert 0.0<=burninfrac<=1.0
    Tacrosschains = gettargtemps(optin)
    iters = TransD_GP.history(optin, stat=:iter)
    start = round(Int, length(iters)*burninfrac)
    start == 0 && (start = 1)
    @info "obtaining models $(iters[start]) to $(iters[end])"
    nmodels = sum((Tacrosschains[start:end,:] .== 1))
    if stat == :fstar || stat == :x_ftrain
        mat1 = Array{Array{Float64}, 1}(undef, nmodels)
    else
        mat1 = Array{Real, 1}(undef, nmodels)
    end
    opt = deepcopy(optin)
    imodel = 0
    for ichain in 1:size(Tacrosschains, 2)
        opt.fstar_filename = "models_"*isns*opt.fdataname*"_$ichain.bin"
        opt.x_ftrain_filename = "points_"*isns*opt.fdataname*"_$ichain.bin"
        opt.costs_filename = "misfits_"*isns*opt.fdataname*"_$ichain.bin"
        if stat == :fstar
            at1idx = findall(Tacrosschains[:,ichain].==1) .>= start
        else
            at1idx = Tacrosschains[:,ichain].==1
            at1idx[1:start-1] .= false
        end
        ninchain = sum(at1idx)
        @info "chain $ichain has $ninchain models"
        ninchain == 0 && continue
        mat1[imodel+1:imodel+ninchain] .= TransD_GP.history(opt, stat=stat)[at1idx]
        imodel += ninchain
    end
    mat1
end

function gettargtemps(opt_in::TransD_GP.Options)
    isns = checkns(opt_in)
    nchains = length(filter( x -> occursin(r"misfits_ns.*bin", x), readdir(pwd()) )) # my terrible regex
    @info "Number of chains is $nchains"
    # now look at any chain to get how many iterations
    opt = deepcopy(opt_in)
    costs_filename = "misfits_"*isns*opt.fdataname
    opt.costs_filename    = costs_filename*"_1.bin"
    iters          = TransD_GP.history(opt, stat=:iter)
    niters         = length(iters)
    @info "McMC has run for $(iters[end]) iterations"
    # then create arrays of unsorted by temperature T
    Tacrosschains  = zeros(Float64, niters, nchains)
    # get the values into the arrays
    for ichain in 1:nchains
        opt_in.costs_filename = costs_filename*"_$ichain.bin"
        Tacrosschains[:,ichain] = TransD_GP.history(opt_in, stat=:T)
    end
    Tacrosschains
end

function checkns(optin::TransD_GP.Options)
    isns = typeof(optin) == TransD_GP.OptionsNonstat
    @info "ns is $isns"
    ns = "ns"
    isns || (ns="s")
    return ns
end

function getchi2forall(opt_in::TransD_GP.Options;
                        nchains          = 1,
                        figsize          = (17,8),
                        fsize            = 14,
                      )
    if nchains == 1 # then actually find out how many chains there are saved
        nchains = length(filter( x -> occursin(r"misfits_ns.*bin", x), readdir(pwd()) )) # my terrible regex
    end
    # now look at any chain to get how many iterations
    isns = checkns(opt_in)
    opt = deepcopy(opt_in)
    costs_filename = "misfits_"*isns*opt.fdataname
    opt_in.costs_filename    = costs_filename*"_1.bin"
    iters          = TransD_GP.history(opt, stat=:iter)
    niters         = length(iters)
    # then create arrays of unsorted by temperature T, k, and chi2
    Tacrosschains  = zeros(Float64, niters, nchains)
    kacrosschains  = zeros(Int, niters, nchains)
    X2by2inchains  = zeros(Float64, niters, nchains)
    # get the values into the arrays
    for ichain in 1:nchains
        opt.costs_filename = costs_filename*"_$ichain.bin"
        Tacrosschains[:,ichain] = TransD_GP.history(opt, stat=:T)
        kacrosschains[:,ichain] = TransD_GP.history(opt, stat=:nodes)
        X2by2inchains[:,ichain] = TransD_GP.history(opt, stat=:U)
    end

    f, ax = plt.subplots(3,2, sharex=true, figsize=figsize)
    ax[1].plot(iters, kacrosschains)
    ax[1].set_title("unsorted by temperature")
    ax[1].grid()
    ax[1].set_ylabel("# nodes")
    ax[2].plot(iters, X2by2inchains)
    ax[2].grid()
    ax[2].set_ylabel("-Log L")
    ax[3].grid()
    ax[3].plot(iters, Tacrosschains)
    ax[3].set_ylabel("Temperature")
    ax[3].set_xlabel("iterations")

    for jstep = 1:niters
        sortidx = sortperm(vec(Tacrosschains[jstep,:]))
        X2by2inchains[jstep,:] = X2by2inchains[jstep,sortidx]
        kacrosschains[jstep,:] = kacrosschains[jstep,sortidx]
        Tacrosschains[jstep,:] = Tacrosschains[jstep,sortidx]
    end

    nchainsatone = sum(Tacrosschains[1,:] .== 1)
    ax[4].plot(iters, kacrosschains)
    ax[4].set_title("sorted by temperature")
    ax[4].plot(iters, kacrosschains[:,1:nchainsatone], "k")
    ax[4].grid()
    ax[5].plot(iters, X2by2inchains)
    ax[5].plot(iters, X2by2inchains[:,1:nchainsatone], "k")
    ax[5].grid()
    ax[6].plot(iters, Tacrosschains)
    ax[6].plot(iters, Tacrosschains[:,1:nchainsatone], "k")
    ax[6].grid()
    ax[6].set_xlabel("iterations")

    nicenup(f, fsize=fsize)

end

function geomprogdepth(n, dy, c)
    dy*(1.0-c^n)/(1-c)
end

function getn(z, dy, c)
    log(1 - z/dy*(1-c))/log(c)
end

function unwrap(v, inplace=false)
  # currently assuming an array
  unwrapped = inplace ? v : copy(v)
  for i in 2:length(v)
    while unwrapped[i] - unwrapped[i-1] >= pi
      unwrapped[i] -= 2pi
    end
    while unwrapped[i] - unwrapped[i-1] <= -pi
      unwrapped[i] += 2pi
    end
  end
  return unwrapped
end

function setupz(zstart, extendfrac;n=100, dz=10)
    @assert extendfrac>1
    znrange        = 1.0:n
    zboundaries    = [zstart, zstart .+ geomprogdepth.(znrange, dz, extendfrac)...]
    thickness      = dz*(extendfrac).^(znrange .-1)
    @assert abs(diff(zboundaries)[end] - thickness[end]) < 1e-12
    zall           = zboundaries[1:end-1] + thickness/2
    znall          = getn.(zall .- zstart, dz, extendfrac)
    plotdepthtransforms(zall, znall, zboundaries, thickness)
    return zall, znall, zboundaries[1:end-1]
end

function setupz(zstart;n=100, dz=10)
    znrange        = 1.0:n
    zboundaries    = range(zstart, step=dz, length=n)
    thickness      = diff(zboundaries)
    zall           = zboundaries[1:end-1] .+ dz/2
    znall          = 1.5:1:n-0.5
    plotdepthtransforms(zall, znall, zboundaries, thickness)
    return zall, znall, zboundaries[1:end-1]
end

function plotdepthtransforms(zall, znall, zboundaries, thickness)
    figure()
    plot(znall, zall)
    xlabel("depth index")
    ylabel("depth associated m")
    grid()
    nicenup(gcf())

    f, ax = plt.subplots(1, 1, figsize=(10,5))
    ax.stem(zboundaries[1:end-1], zboundaries[1:end-1], markerfmt="")
    ax.stem(zall, zall, "k--", markerfmt=" ")
    ax.set_xlabel("depth m")
    ax.set_ylabel("depth m")
    nicenup(gcf())

    f, ax = plt.subplots(1, 2, sharey=true, figsize=(11,5))
    ax[1].stem(zall,thickness, "k--", markerfmt=" ")
    ax[1].set_ylabel("thickness m")
    ax[2].set_xlabel("depth m")
    ax[1].yaxis.grid(which="major")
    ax[2].stem(znall,thickness, "k--", markerfmt=" ")
    ax[2].set_xlabel("depth index")
    ax[2].yaxis.grid(which="major")
    nicenup(f)
end

function nicenup(g::PyPlot.Figure;fsize=16)
    for ax in gcf().axes
        ax.tick_params("both",labelsize=fsize)
        ax.xaxis.label.set_fontsize(fsize)
        ax.yaxis.label.set_fontsize(fsize)
        any(keys(ax) .== :zaxis) && ax.zaxis.label.set_fontsize(fsize)
        ax.title.set_fontsize(fsize)

        if typeof(ax.get_legend_handles_labels()[1]) != Array{Any,1}
            ax.legend(loc="best", fontsize=fsize)
        end
    end
    g.tight_layout()
end

function plot_posterior(operator::Operator1D,
                        optns::TransD_GP.OptionsNonstat,
                        opt::TransD_GP.OptionsStat;
    nbins = 50,
    burninfrac=0.5,
    isns=true,
    qp1=0.05,
    qp2=0.95,
    cmappdf = "viridis",
    figsize=(10,5),
    pdfnormalize=false,
    fsize=14)
    himage_ns, edges_ns, CI_ns = make1Dhist(optns, burninfrac=burninfrac, nbins = nbins, qp1=qp1, qp2=qp2, pdfnormalize=pdfnormalize)
    himage, edges, CI = make1Dhist(opt, burninfrac=burninfrac, nbins = nbins, qp1=qp1, qp2=qp2, pdfnormalize=pdfnormalize)
    f,ax = plt.subplots(1,2, sharey=true, figsize=figsize)
    xall = opt.xall
    im1 = ax[1].imshow(himage_ns, extent=[edges_ns[1],edges_ns[end],xall[end],xall[1]], aspect="auto", cmap=cmappdf)
    cb1 = colorbar(im1, ax=ax[1])
    cb1.ax.set_xlabel("pdf \nns")
    ax[1].grid()
    im2 = ax[2].imshow(himage, extent=[edges[1],edges[end],xall[end],xall[1]], aspect="auto", cmap=cmappdf)
    ax[2].grid()
    cb2 = colorbar(im2, ax=ax[2])
    cb2.ax.set_xlabel("pdf \nstationary")
    ax[1].plot(CI_ns, xall[:], linewidth=2, color="g")
    ax[1].plot(CI_ns, xall[:], linewidth=2, color="k", linestyle="--")
    ax[2].plot(CI, xall[:], linewidth=2, color="g")
    ax[2].plot(CI, xall[:], linewidth=2, color="k", linestyle="--")
    ax[1].set_xlabel(L"\log_{10} \rho")
    ax[1].set_ylabel("depth (m)")
    ax[2].set_xlabel(L"\log_{10} λ")
    nicenup(f, fsize=fsize)
end

function plot_posterior(operator::Operator1D,
                        opt::TransD_GP.OptionsStat;
    nbins = 50,
    burninfrac=0.5,
    isns=true,
    qp1=0.05,
    qp2=0.95,
    cmappdf = "viridis",
    figsize=(5,5),
    pdfnormalize=false,
    fsize=14)
    himage, edges, CI = make1Dhist(opt, burninfrac=burninfrac, nbins = nbins, qp1=qp1, qp2=qp2, pdfnormalize=pdfnormalize)
    f,ax = plt.subplots(1,1, sharey=true, figsize=figsize)
    xall = opt.xall
    #im1 = ax.imshow(himage, extent=[edges[1],edges[end],xall[end],xall[1]], aspect="auto", cmap=cmappdf)
    #ax.grid()
    im1 = ax.pcolormesh(edges[:], xall[:], himage, cmap=cmappdf)
    cb1 = colorbar(im1, ax=ax)
    cb1.ax.set_xlabel("pdf \nstationary")
    ax.plot(CI, xall[:], linewidth=2, color="g")
    ax.plot(CI, xall[:], linewidth=2, color="k", linestyle="--")
    ax.set_xlabel(L"\log_{10} \rho")
    ax.set_ylabel("depth (m)")
    nicenup(f, fsize=fsize)
end

function make1Dhist(opt::TransD_GP.Options;
                burninfrac = 0.5,
                nbins = 50,
                rhomin=Inf,
                rhomax=-Inf,
                qp1=0.05,
                qp2=0.95,
                pdfnormalize=false)
    M = assembleTat1(opt, :fstar, burninfrac=burninfrac)
    if (rhomin == Inf) && (rhomax == -Inf)
        for (i,mm) in enumerate(M)
            rhomin_mm = minimum(mm)
            rhomax_mm = maximum(mm)
            rhomin_mm < rhomin && (rhomin = rhomin_mm)
            rhomax_mm > rhomax && (rhomax = rhomax_mm)

        end
        if typeof(opt) == TransD_GP.OptionsStat && opt.needλ²fromlog
            rhomin = 0.5*log10(rhomin)
            rhomax = 0.5*log10(rhomax)
        end
    else
        @assert rhomin < rhomin
    end
    edges = LinRange(rhomin, rhomax, nbins+1)
    himage = zeros(Float64, length(M[1]), nbins)
    CI = zeros(Float64, length(M[1]), 2)
    for ilayer=1:size(opt.xall,2)
        if typeof(opt) == TransD_GP.OptionsStat && opt.needλ²fromlog
            himage[ilayer,:] = fit(Histogram, [0.5log10.(m[ilayer]) for m in M], edges).weights
        else
            himage[ilayer,:] = fit(Histogram, [m[ilayer] for m in M], edges).weights
        end
        himage[ilayer,:] = himage[ilayer,:]/sum(himage[ilayer,:])/(diff(edges)[1])
        pdfnormalize && (himage[ilayer,:] = himage[ilayer,:]/maximum(himage[ilayer,:]))
        if typeof(opt) == TransD_GP.OptionsStat && opt.needλ²fromlog
            CI[ilayer,:] = [quantile([0.5log10.(m[ilayer]) for m in M],(qp1, qp2))...]
        else
            CI[ilayer,:] = [quantile([m[ilayer] for m in M],(qp1, qp2))...]
        end
    end
    return himage, edges, CI
end

function make1Dhists(opt_in::TransD_GP.Options, burninfrac::Real;
                        kfigsize=(8,8),
                        nxbins=50,
                        nftbins=50,
                        qp1=0.01,
                        qp2=0.99)
    f2, ax2 = plt.subplots(2,2, figsize=kfigsize)
    x, ft = trimxft(opt_in, burninfrac)
    n = assembleTat1(opt_in, :nodes, burninfrac=burninfrac)
    edgesx = LinRange(opt_in.xbounds[1], opt_in.xbounds[2], nxbins+1)
    edgesrho = LinRange(opt_in.fbounds[1], opt_in.fbounds[2], nftbins+1)
    h = fit(Histogram, (x[:], ft[:]), (edgesx, edgesrho)).weights
    vmin, vmax = quantile(h[:], (qp1, qp2))
    im3 = ax2[1].imshow(h, extent=[edgesrho[1],edgesrho[end],edgesx[end],edgesx[1]], aspect="auto", vmin=vmin, vmax=vmax)
    ax2[1].set_ylabel("depth m")
    ax2[1].set_xlabel(L"\log_{10}\rho")
    rhist = sum(h,dims=1)[:]
    rhist = rhist/sum(rhist)/diff(edgesrho)[1]
    ax2[2].bar(0.5*(edgesrho[2:end]+edgesrho[1:end-1]), rhist, width=diff(edgesrho)[1])
    ax2[2].plot([edgesrho[1], edgesrho[end]], 1/(opt_in.fbounds[end]-opt_in.fbounds[1])*[1,1],"--k")
    ax2[2].set_ylim(0, 3*maximum(rhist))
    ax2[2].set_xlabel(L"\log_{10}\rho")
    ax2[2].set_ylabel("pdf")
    xhist = sum(h,dims=2)[:]
    xhist = xhist/sum(xhist)/diff(edgesx)[1]
    ax2[3].barh(0.5*(edgesx[2:end]+edgesx[1:end-1]), xhist, height=diff(edgesx)[1])
    ax2[3].plot(1/(opt_in.xbounds[end]-opt_in.xbounds[1])*[1,1], [edgesx[1], edgesx[end]],"--k")
    ax2[3].set_xlim(0, 3*maximum(xhist))
    ax2[3].set_ylim(ax2[1].get_ylim())
    ax2[3].set_ylabel("depth m")
    ax2[3].set_xlabel("pdf")
    # ax2[3][:xaxis][:set_label_position]("top")
    ax2[2].get_shared_y_axes().join(ax2[1], ax2[3])
    ax2[2].get_shared_x_axes().join(ax2[1], ax2[2])
    ax2[2].set_xlim(ax2[1].get_xlim())
    k = fit(Histogram, n, (opt_in.nmin-0.5:opt_in.nmax+0.5)).weights
    k = k/sum(k)
    @show mean(k)
    ax2[4].bar(opt_in.nmin:opt_in.nmax, k, width=1)
    ax2[4].plot([opt_in.nmin, opt_in.nmax],1/(opt_in.nmax-opt_in.nmin+1)*[1,1],"--k")
    ax2[4].set_xlabel("# training points")
    ax2[4].set_ylim(0, 3*max(k...))
    ax2[4].set_ylabel("pdf")
    ax2[4].set_xlim(opt_in.nmin, opt_in.nmax)
end

end
