module CommonToAll
using PyPlot, StatsBase, Statistics, Distances, LinearAlgebra,
      DelimitedFiles, ..AbstractOperator, NearestNeighbors, CSV, DataFrames

import ..Options, ..OptionsStat, ..OptionsNonstat, ..OptionsNuisance,
       ..history, ..GP.κ, ..calcfstar!, ..AbstractOperator.Sounding

export trimxft, assembleTat1, gettargtemps, checkns, getchi2forall, nicenup,
        plot_posterior, make1Dhist, make1Dhists, setupz, zcontinue, makezρ, plotdepthtransforms,
        unwrap, getn, geomprogdepth, assemblemodelsatT, getstats, gethimage,
        assemblenuisancesatT, makenuisancehists, stretchexists,
        makegrid, whichislast, makesummarygrid, makearray, plotNEWSlabels, 
        plotprofile, gridpoints, splitsoundingsbyline, dfn2hdr

function trimxft(opt::Options, burninfrac::Float64, temperaturenum::Int)
    x_ft = assembleTat1(opt, :x_ftrain, burninfrac=burninfrac, temperaturenum=temperaturenum)
    n = assembleTat1(opt, :nodes, burninfrac=burninfrac, temperaturenum=temperaturenum)
    x, ft = zeros(size(opt.xall, 1), sum(n)), zeros(size(opt.fbounds, 1), sum(n))
    nlast = 0
    for (i, xft) in enumerate(x_ft)
        x[:,nlast+1:nlast+n[i]]  = xft[1:size(opt.xall, 1), 1:n[i]]
        ft[:,nlast+1:nlast+n[i]] = xft[size(opt.xall, 1)+1:end, 1:n[i]]
        nlast += n[i]
    end
    x, ft, n
end

function assembleTat1(optin::Options, stat::Symbol; burninfrac=0.5, temperaturenum=-1)
    @assert temperaturenum!=-1 "Please explicitly specify which temperature number"
    if stat == :fstar && temperaturenum!= 1
        return assemblemodelsatT(optin, burninfrac=burninfrac, temperaturenum=temperaturenum)
    end
    isns = checkns(optin)
    @assert 0.0<=burninfrac<=1.0
    Tacrosschains = gettargtemps(optin)
    temps = sort(unique(Tacrosschains))
    niters = size(Tacrosschains, 1)
    start = round(Int, niters*burninfrac)
    start == 0 && (start = 1)
    ttarget = temps[temperaturenum]
    nmodels = sum((Tacrosschains[start:end,:] .== ttarget))
    if stat == :fstar || stat == :x_ftrain
        mat1 = Array{Array{Float64}, 1}(undef, nmodels)
    else
        mat1 = Array{Real, 1}(undef, nmodels)
    end
    opt = deepcopy(optin)
    imodel = 0
    for ichain in 1:size(Tacrosschains, 2)
        opt.fstar_filename = "models_"*opt.fdataname*isns*"_$ichain.bin"
        opt.x_ftrain_filename = "points_"*opt.fdataname*isns*"_$ichain.bin"
        opt.costs_filename = "misfits_"*opt.fdataname*isns*"_$ichain.bin"
        if stat == :fstar
            at1idx = findall(Tacrosschains[:,ichain].==ttarget) .>= start
        else
            at1idx = Tacrosschains[:,ichain].==ttarget
            at1idx[1:start-1] .= false
        end
        ninchain = sum(at1idx)
        @info "chain $ichain has $ninchain models"
        ninchain == 0 && continue
        mat1[imodel+1:imodel+ninchain] .= history(opt, stat=stat)[at1idx]
        imodel += ninchain
    end
    iters = history(opt, stat=:iter)
    @info "obtained models $(iters[start]) to $(iters[end]) at T=$ttarget"
    mat1
end

function assemblemodelsatT(opt::OptionsStat; burninfrac=0.9, temperaturenum=-1)
    @assert temperaturenum!=-1 "Please explicitly specify which temperature number"
    x_ft = assembleTat1(opt, :x_ftrain, burninfrac=burninfrac, temperaturenum=temperaturenum)
    n = assembleTat1(opt, :nodes, burninfrac=burninfrac, temperaturenum=temperaturenum)
    nmodels = length(n)
    matT = Array{Array{Float64}, 1}(undef, nmodels)
    K_y = zeros(opt.nmax, opt.nmax)
    Kstar = zeros(Float64, size(opt.xall,2), opt.nmax)
    xtest = opt.xall
    for imodel = 1:nmodels
        xtrain = @view x_ft[imodel][1:size(opt.xall, 1), 1:n[imodel]]
        ftrain = x_ft[imodel][size(opt.xall, 1)+1:end, :]
        ky = view(K_y, 1:n[imodel], 1:n[imodel])
        map!(x->κ(opt.K, x),ky,pairwise(WeightedEuclidean(1 ./opt.λ² ), xtrain, dims=2))
        K_y[diagind(K_y)] .+= opt.δ^2
        ks = view(Kstar, :, 1:n[imodel])
        map!(x->κ(opt.K, x),Kstar,pairwise(WeightedEuclidean(1 ./opt.λ² ), xtest, xtrain, dims=2))
        matT[imodel] = zeros(size(opt.fbounds, 1), size(opt.xall, 2))
        fstar = matT[imodel]
        calcfstar!(fstar, ftrain, opt, K_y, Kstar, n[imodel])
    end
    matT
end

function assemblenuisancesatT(optn::OptionsNuisance;
    burninfrac = 0.5, temperaturenum = -1)
    @assert temperaturenum != -1 "Please specify temperature idx explicitly"
    @assert 0.0 <= burninfrac < 1.0
    Tacrosschains = gettargtemps(optn)
    #this is probably insanely inefficient
    #Θ(niters*nchains) to run the unique() call.
    #however, we need to iterate over the entire array
    #at some point to build the ensemble so it's
    #fine. Possible improvement: store temperatures
    #as the "temperature number" to start with and keep
    #a separate mapping from temperature index to the float
    #value. This means our sortedTs is just 1:nchains and
    #also avoids any issues with equality testing of floats.
    sortedTs = sort(unique(Tacrosschains))
    niters = size(Tacrosschains,1)
    #this will never give a bounds error
    #because of the assert above
    firsti = round(Int, niters*burninfrac)
    firsti == 0 && (start = 1)
    # firsti = 1 + floor(Int, niters*burninfrac)
    ttarg = sortedTs[temperaturenum]
    nmodels = sum(Tacrosschains[firsti:end,:] .== ttarg)
    fdataname = optn.fdataname
    #drop iteration number
    mvals = zeros(nmodels,optn.nnu)
    imodel = 0
    for ichain = 1:size(Tacrosschains,2)
        at1idx = Tacrosschains[:,ichain].==ttarg
        at1idx[1:firsti-1] .= false
        ninchain = sum(at1idx)
        @info "chain $ichain has $ninchain models"
        ninchain == 0 && continue
        vals_filename = "values_nuisance_"*fdataname*"$ichain.bin"
        ndat = readdlm(vals_filename, ' ', Float64)[:,2:end]
        mvals[imodel+1:imodel+ninchain,:] .= ndat[at1idx,:]
        imodel += ninchain
    end
    mvals
end

function getnchains(costs_filename)
c = 0
    for fname in readdir(pwd())
           r = Regex(costs_filename)
           c += !isa(match(r, fname), Nothing)
    end
    c
end

function gettargtemps(opt_in::Options)
    isns = checkns(opt_in)
    opt = deepcopy(opt_in)
    costs_filename = "misfits_"*opt.fdataname*isns
    nchains = getnchains(costs_filename)
    @info "Number of chains is $nchains"
    # now look at any chain to get how many iterations
    opt.costs_filename    = costs_filename*"_1.bin"
    iters          = history(opt, stat=:iter)
    niters         = length(iters)
    @info "McMC has run for $(iters[end]) iterations"
    # then create arrays of unsorted by temperature T
    Tacrosschains  = zeros(Float64, niters, nchains)
    # get the values into the arrays
    for ichain in 1:nchains
        opt.costs_filename = costs_filename*"_$ichain.bin"
        Tacrosschains[:,ichain] = history(opt, stat=:T)
    end
    Tacrosschains
end

function gettargtemps(optn_in::OptionsNuisance)
    optn = deepcopy(optn_in)
    costs_filename = "misfits_"*optn.fdataname*"nuisance"
    nchains = getnchains(costs_filename)
    @info "Number of chains is $nchains"
    optn.costs_filename = costs_filename*"_1.bin"
    iters = history(optn, stat=:iter)
    niters = length(iters)
    @info "MCMC has run for $(iters[end]) iterations"
    Tacrosschains = zeros(Float64, niters, nchains)
    for ichain in 1:nchains
        optn.costs_filename = costs_filename*"_$ichain.bin"
        Tacrosschains[:,ichain] = history(optn, stat=:T)
    end
    Tacrosschains
end

function getstats(optin::Options;
                  figsize=(5,6), fontsize=12,
                  nxticks=5, nyticks=5, alpha=0.6)
    isns = checkns(optin)
    opt = deepcopy(optin)
    costs_filename = "misfits_"*opt.fdataname*isns
    nchains = getnchains(costs_filename)
    chains = 1:nchains
    @info "Number of chains is $nchains"
    opt.costs_filename    = costs_filename*"_1.bin"
    iters          = history(opt, stat=:iter)
    statnames = [:acceptanceRateBirth, :acceptanceRateDeath,
                 :acceptanceRatePosition, :acceptanceRateProperty, :acceptanceRateDC]
    f,ax = plt.subplots(length(statnames), 1,
                        sharex=true, sharey=true, figsize=figsize)
    maxar = 0
    for (ichain, chain) in enumerate(chains)
        opt.costs_filename = costs_filename*"_$chain.bin"
        for (istat, statname) in enumerate(statnames)
            ar = history(opt, stat=statname)
            mx = maximum(ar[.!isnan.(ar)])
            mx > maxar && (maxar = mx)
            ax[istat].plot(iters, ar, alpha=alpha)
            if ichain == nchains
                ax[istat].set_title("$statname "*isns)
                ax[istat].set_ylabel("acc. %")
                if istat == length(statnames)
                    ax[istat].set_xlabel("mcmc step no.")
                    ax[istat].set_xticks(LinRange(iters[1], iters[end], nxticks))
                    ax[istat].set_xlim(extrema(iters))
                    ax[istat].set_yticks(LinRange(0, maxar, nyticks))
                    ax[istat].set_ylim(0, maxar)
                end
            end
        end
    end
    nicenup(f, fsize=fontsize)
end

function getstats(optin::OptionsNuisance;
                  figsize=(5,6), fontsize=12,
                  nxticks=5, nyticks=5, alpha=0.6)
    opt = deepcopy(optin)
    costs_filename = "misfits_"*opt.fdataname*"nuisance"
    nchains = getnchains(costs_filename)
    chains = 1:nchains
    @info "Number of chains is $nchains"
    opt.costs_filename = costs_filename*"_1.bin"
    iters          = history(opt, stat=:iter)
    statname = :acceptanceRate
    f,ax = plt.subplots(length(optin.idxnotzero), 1, sharex=true, sharey=true, figsize=figsize)
    maxar = 0
    for (ichain, chain) in enumerate(chains)
        opt.costs_filename = costs_filename*"_$chain.bin"
        ar = history(opt, stat=statname)[:,optin.idxnotzero]
        for (i, idx) in enumerate(optin.idxnotzero)
            mx = maximum(ar[.!isnan.(ar)])
            mx > maxar && (maxar = mx)
            ax[i].plot(iters, ar[:,i], alpha=alpha)
            if ichain == nchains
                ax[i].set_title("Nuisance #$idx")
                ax[i].set_ylabel("acc. %")
                if i == length(optin.idxnotzero)
                    ax[i].set_xlabel("mcmc step no.")
                    ax[i].set_xticks(LinRange(iters[1], iters[end], nxticks))
                    ax[i].set_xlim(extrema(iters))
                end
            end
        end
    end
    nicenup(f, fsize=fontsize)
end

function checkns(optin::Options)
    isns = typeof(optin) == OptionsNonstat
    @info "ns is $isns"
    ns = "ns"
    isns || (ns="s")
    return ns
end

function getchi2forall(optn_in::OptionsNuisance;
                       nchains = 1,
                       figsize = (6,4),
                       fsize = 8,
                       alpha = 0.25,
                       nxticks = 0,
                       gridon = false)
    optn = deepcopy(optn_in)
    fdataname = optn.fdataname
    costs_filename = "misfits_"*fdataname*"nuisance"
    if nchains == 1
        nchains = getnchains(costs_filename)
    end
    optn.costs_filename = costs_filename*"_1.bin"
    iters = history(optn, stat=:iter)
    niters = length(iters)
    Tacrosschains = zeros(Float64, niters, nchains)
    χ2acrosschains = zeros(Float64, niters, nchains)
    for chain in 1:nchains
        optn.costs_filename = costs_filename*"_$(chain).bin"
        Tacrosschains[:,chain] = history(optn, stat=:T)
        χ2acrosschains[:,chain] = history(optn, stat=:misfit)
    end
    Torder = sort([(i,j) for i=1:niters, j=1:nchains],
                    by = ix->Tacrosschains[ix...],
                    dims = 2)
    χ2sorted = [χ2acrosschains[Torder[i,j]...] for i=1:niters, j=1:nchains]
    fig, ax = subplots(3,1, sharex=true, figsize=figsize)
    ax[1].plot(iters, χ2acrosschains, alpha=alpha)
    ax[1].set_xlim(extrema(iters)...)
    ax[1].set_ylim(0,100)
    ax[1].set_title("unsorted χ^2 misfit")
    ax[2].plot(iters, Tacrosschains, alpha=alpha)
    ax[2].set_title("temperature")
    ax[3].plot(iters, χ2sorted, alpha=alpha)
    ax[3].set_ylim(0,100)
    ax[3].set_title("sorted χ^2 misfit")
    ax[3].set_xlabel("iterations")

    nxticks == 0 || ax[3].set_xticks(iters[1]:div(iters[end],nxticks):iters[end])
    nicenup(fig, fsize = fsize)
end

function getchi2forall(opt_in::Options;
                        nchains         = 1,
                        figsize         = (6,4),
                        fsize           = 8,
                        alpha           = 0.25,
                        nxticks         = 0,
                        gridon          = false,
                        omittemp        = false,
                      )
    # now look at any chain to get how many iterations
    isns = checkns(opt_in)
    opt = deepcopy(opt_in)
    costs_filename = "misfits_"*opt.fdataname*isns
    if nchains == 1 # then actually find out how many chains there are saved
        nchains = getnchains(costs_filename)
    end
    opt.costs_filename    = costs_filename*"_1.bin"
    iters          = history(opt, stat=:iter)
    niters         = length(iters)
    # then create arrays of unsorted by temperature T, k, and chi2
    Tacrosschains  = zeros(Float64, niters, nchains)
    kacrosschains  = zeros(Int, niters, nchains)
    X2by2inchains  = zeros(Float64, niters, nchains)
    # get the values into the arrays
    for ichain in 1:nchains
        opt.costs_filename = costs_filename*"_$ichain.bin"
        Tacrosschains[:,ichain] = history(opt, stat=:T)
        kacrosschains[:,ichain] = history(opt, stat=:nodes)
        X2by2inchains[:,ichain] = history(opt, stat=:U)
    end

    f, ax = plt.subplots(3,1, sharex=true, figsize=figsize)
    ax[1].plot(iters, kacrosschains, alpha=alpha)
    ax[1].set_xlim(extrema(iters)...)
    ax[1].set_title(isns*" unsorted by temperature")
    gridon && ax[1].grid()
    ax[1].set_ylabel("# nodes")
    ax[2].plot(iters, X2by2inchains, alpha=alpha)
    gridon && ax[2].grid()
    ax[2].set_ylabel("-Log L")
    gridon && ax[3].grid()
    ax[3].plot(iters, Tacrosschains, alpha=alpha)
    ax[3].set_ylabel("Temperature")
    ax[3].set_xlabel("iterations")
    nxticks == 0 || ax[3].set_xticks(iters[1]:div(iters[end],nxticks):iters[end])
    nicenup(f, fsize=fsize)

    Tunsorted = copy(Tacrosschains)
    for jstep = 1:niters
        sortidx = sortperm(vec(Tacrosschains[jstep,:]))
        X2by2inchains[jstep,:] = X2by2inchains[jstep,sortidx]
        kacrosschains[jstep,:] = kacrosschains[jstep,sortidx]
        Tacrosschains[jstep,:] = Tacrosschains[jstep,sortidx]
    end
    nrows = omittemp ? 2 : 3
    f, ax = plt.subplots(nrows, 1, sharex=true, figsize=figsize)
    nchainsatone = sum(Tacrosschains[1,:] .== 1)
    ax[1].plot(iters, kacrosschains, alpha=alpha)
    ax[1].set_xlim(extrema(iters)...)
    ax[1].set_ylabel("# nuclei")
    ax[1].set_title(isns*" sorted by temperature")
    ax[1].plot(iters, kacrosschains[:,1:nchainsatone], "k", alpha=alpha)
    gridon && ax[1].grid()
    ax[2].plot(iters, X2by2inchains, alpha=alpha)
    ax[2].plot(iters, X2by2inchains[:,1:nchainsatone], "k", alpha=alpha)
    ax[2].set_ylabel("-Log L")
    gridon && ax[2].grid()
    if !omittemp
        ax[3].plot(iters, Tunsorted, alpha=alpha, color="gray")
        ax[3].set_ylabel("Temperature")
        gridon && ax[3].grid()
    end
    nxticks == 0 || ax[3].set_xticks(iters[1]:div(iters[end],nxticks):iters[end])
    ax[nrows].set_xlabel("iterations")
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

function setupz(zstart, extendfrac;n=100, dz=10, showplot=false, forextension=false)
    @assert extendfrac>1
    znrange        = 1.0:n
    zboundaries    = [zstart, zstart .+ geomprogdepth.(znrange, dz, extendfrac)...]
    thickness      = diff(zboundaries)
    zall           = zboundaries[1:end-1] + thickness/2
    znall          = getn.(zall .- zstart, dz, extendfrac)
    showplot && plotdepthtransforms(zall, znall, zboundaries)
    zbreturn = forextension ? zboundaries : zboundaries[1:end-1]
    return zall, znall, zbreturn
end

function zcontinue(;zall=nothing, znall=nothing, zboundaries=nothing, n=nothing,
                extendfrac=nothing, dz=nothing, showplot=false)
    isa(dz, Nothing) && (dz=diff(zboundaries)[end])
    zall_, znall_, zboundaries_ = setupz(zboundaries[end], extendfrac, n=n, dz=dz, showplot=false)
    zall, znall, zboundaries = [zall;zall_], [znall; znall[end] .+ znall_], [zboundaries[1:end-1]; zboundaries_]
    zlast = zboundaries[end] + diff(zboundaries)[end]
    showplot && plotdepthtransforms(zall, znall, [zboundaries;zlast])
    @assert all(diff(zall).>0)
    @assert all(diff(znall).>0)
    @assert all(diff(zboundaries).>0)
    zall, znall, zboundaries
end

function setupz(zstart;n=100, dz=10)
    zboundaries    = range(zstart, step=dz, length=n)
    zall           = zboundaries[1:end-1] .+ dz/2
    znall          = 1.5:1:n-0.5
    plotdepthtransforms(zall, znall, zboundaries)
    return zall, znall, zboundaries[1:end-1]
end

function plotdepthtransforms(zall, znall, zboundaries)
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

    thickness = diff(zboundaries)
    f, ax = plt.subplots(1, 2, sharey=true, figsize=(11,5))
    ax[1].stem(zall,thickness, "k--", markerfmt=" ")
    ax[1].set_ylabel("thickness m")
    ax[1].set_xlabel("depth m")
    ax[1].yaxis.grid(which="major")
    ax[2].stem(znall,thickness, "k--", markerfmt=" ")
    ax[2].set_xlabel("depth index")
    ax[2].yaxis.grid(which="major")
    nicenup(gcf())
end

function makezρ(zboundaries::Array{Float64, 1};
          zfixed      = [-1e6,    0.],
          ρfixed      = [1e12     0.3])
    @assert length(zfixed) == length(ρfixed)
    z = [zfixed..., zboundaries...]
    @assert all(diff(z).>0)
    ρ = zeros(length(z))
    nfixed = length(ρfixed)
    ρ[1:nfixed] .= ρfixed
    z, ρ, nfixed
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

function plot_posterior(F::Operator1D,
                        optns::OptionsNonstat,
                        opts::OptionsStat;
                        temperaturenum = 1,
                        nbins = 50,
                        burninfrac=0.5,
                        qp1=0.05,
                        qp2=0.95,
                        vmaxpc=1.0,
                        cmappdf = "inferno",
                        figsize=(12,5),
                        pdfnormalize=false,
                        fsize=14,
                        showlscale1vd=false,
                        istothepow=false,
                        doplot = true,
                        a = 0.25,
                        lw = 1)
    @assert 0<vmaxpc<=1
    M = assembleTat1(optns, :fstar, burninfrac=burninfrac, temperaturenum=temperaturenum)
    himage_ns, edges_ns, CI_ns, meanimage_ns, 
    meandiffimage_ns, sdslope_ns, = gethimage(F, M, optns, temperaturenum=temperaturenum,
                            nbins=nbins, qp1=qp1, qp2=qp2, istothepow=istothepow,
                            islscale=false, pdfnormalize=pdfnormalize)
    M = assembleTat1(opts, :fstar, burninfrac=burninfrac, temperaturenum=temperaturenum)
    himage, edges, CI, meanimage, 
    meandiffimage, sdslope, = gethimage(F, M, opts, temperaturenum=temperaturenum,
                            nbins=nbins, qp1=qp1, qp2=qp2, istothepow=false,
                            islscale=true, pdfnormalize=pdfnormalize)
    f,ax = plt.subplots(1, 3+convert(Int, showlscale1vd), sharey=true, figsize=figsize)
    xall = opts.xall
    diffs = diff(xall[:])
    xmesh = vcat(xall[1:end-1] - diffs/2, xall[end]-diffs[end]/2, xall[end])
    vmin, vmax = extrema(himage_ns)
    vmax = vmin+vmaxpc*(vmax-vmin)
    im1 = ax[1].pcolormesh(edges_ns[:], xmesh, himage_ns, cmap=cmappdf, vmax=vmax)
    cb1 = colorbar(im1, ax=ax[1])
    cb1.ax.set_xlabel("pdf \nns")
    @show propmin, propmax = getbounds(CI_ns, optns.fbounds)
    ax[1].set_xlim(propmin, propmax)
    ax[2].fill_betweenx(xall[:],meandiffimage_ns[:]-sdslope_ns[:],meandiffimage_ns[:]+sdslope_ns[:],alpha=a, color="gray")
    ax[2].plot(meandiffimage_ns, xall[:], linewidth=lw, color="k")
    ax[2].set_xlabel("slope")

    vmin, vmax = extrema(himage)
    vmax = vmin+vmaxpc*(vmax-vmin)
    im3 = ax[3].pcolormesh(edges[:], xmesh, himage, cmap=cmappdf, vmax=vmax)
    ax[3].set_ylim(extrema(xall)...)
    propmin, propmax = getbounds(CI, opts.fbounds)
    ax[3].set_xlim(propmin, propmax)
    ax[3].invert_yaxis()
    cb3 = colorbar(im3, ax=ax[3])
    cb3.ax.set_xlabel("pdf \nstationary")

    if showlscale1vd
        ax[4].plot(meandiffimage[:], xall[:], linewidth=2, color="g", alpha=0.5)
        ax[4].fill_betweenx(xall[:], meandiffimage[:]-sdslope[:],meandiffimage[:]+sdslope[:], alpha=.5)
        ax[4].set_xlabel("slope")
    end

    ax[1].plot(CI_ns, xall[:], linewidth=lw, color="w", alpha=a)
    # ax[1].plot(CI_ns, xall[:], linewidth=1, color="c", linestyle="--", alpha=0.5)
    # ax[1].plot(meanimage_ns[:], xall[:], linewidth=2, color="w", alpha=0.5)
    # ax[1].plot(meanimage_ns[:], xall[:], linewidth=2, color="k", linestyle="--", alpha=0.5)
    ax[3].plot(CI, xall[:], linewidth=lw, color="w", alpha=a)
    # ax[3].plot(CI, xall[:], linewidth=1, color="c", linestyle="--", alpha=0.5)
    # ax[3].plot(meanimage[:], xall[:], linewidth=2, color="w", alpha=0.5)
    # ax[3].plot(meanimage[:], xall[:], linewidth=2, color="k", linestyle="--", alpha=0.5)
    ax[1].set_xlabel(L"\log_{10} \rho")
    ax[1].set_ylabel("depth (m)")
    ax[3].set_xlabel(L"\log_{10} λ")
    nicenup(f, fsize=fsize)
end

function plot_posterior(F::Operator1D,
                        opt::OptionsStat;
                    temperaturenum = 1,
                    nbins = 50,
                    burninfrac=0.5,
                    qp1=0.05,
                    qp2=0.95,
                    vmaxpc=1.0,
                    cmappdf = "inferno",
                    figsize=(5,5),
                    pdfnormalize=false,
                    istothepow=false,
                    fsize=14,
                    plotCI = true,
                    plotmean = true,
                    alpha = 1.0,
                    CIcolor = ["w", "k"],
                    meancolor = ["m", "r"],
                    lwidth = 2,
                    doplot = true)
    @assert 0<vmaxpc<=1
    
    M = assembleTat1(opt, :fstar, burninfrac=burninfrac, temperaturenum=temperaturenum)
    himage, edges, CI, meanimage, meandiffimage, sdslope, = gethimage(F, M, opt, burninfrac=burninfrac, temperaturenum=temperaturenum,
                nbins=nbins, qp1=qp1, qp2=qp2, istothepow=istothepow,
                islscale=false, pdfnormalize=pdfnormalize)

    if doplot
        f,ax = plt.subplots(1,2, sharey=true, figsize=figsize)
        xall = opt.xall
        diffs = diff(xall[:])
        xmesh = vcat(xall[1:end-1] - diffs/2, xall[end]-diffs[end]/2, xall[end])
        vmin, vmax = extrema(himage)
        vmax = vmin+vmaxpc*(vmax-vmin)
        im1 = ax[1].pcolormesh(edges[:], xmesh, himage, cmap=cmappdf, vmax=vmax)
        ax[1].set_ylim(extrema(xall)...)
        ax[1].invert_yaxis()
        cb1 = colorbar(im1, ax=ax[1])
        cb1.ax.set_xlabel("pdf \nstationary")
        plotCI && ax[1].plot(CI, xall[:], linewidth=lwidth, color=CIcolor[1], alpha=alpha)
        plotCI && ax[1].plot(CI, xall[:], linewidth=lwidth, color=CIcolor[2], linestyle="--", alpha=alpha)
        plotmean && ax[1].plot(meanimage[:], xall[:], linewidth=lwidth, color=meancolor[1], alpha=alpha)
        plotmean && ax[1].plot(meanimage[:], xall[:], linewidth=lwidth, color=meancolor[2], linestyle="--", alpha=alpha)
        ax[1].set_xlabel(L"\log_{10} \rho")
        ax[1].set_ylabel("depth (m)")
        bounds = copy(opt.fbounds)
        if stretchexists(F)
            bounds = [minimum(F.low) maximum(F.low + F.Δ)]
        end
        istothepow && (bounds .= 10 .^ bounds)
        propmin, propmax = getbounds(CI, bounds)
        ax[1].set_xlim(propmin, propmax)
        ax[2].plot(meandiffimage[:], xall[:], linewidth=2, color="k", linestyle="-")
        zeroside = meandiffimage[:]-sdslope[:]
        zeroside[zeroside .< 0] .= 0
        ax[2].fill_betweenx(xall[:],zeroside,meandiffimage[:]+sdslope[:],alpha=.25)
        ax[2].set_xlabel("mean slope")
        nicenup(f, fsize=fsize)
    end
    CI[:,1], CI[:,2], CI[:,3], meanimage, meandiffimage, sdslope
end

function getbounds(CI, bounds)
    propmin = min(minimum(CI), minimum(bounds))
    propmax = max(maximum(CI), maximum(bounds))
    propmin, propmax
end    

function plot_posterior(operator::Operator1D,
                        optn::OptionsNuisance;
                        temperaturenum = 1,
                        nbins = 50,
                        burninfrac=0.5,
                        figsize=(8,16),
                        pdfnormalize=false,
                        fsize=14,
                        qp1 = 0.05,
                        qp2 = 0.95,
                        doplot=true)
    hists, CI = makenuisancehists(optn, qp1, qp2, burninfrac = burninfrac,
                                     nbins = nbins, temperaturenum = temperaturenum,
                                     pdfnormalize = pdfnormalize)
    if doplot
        fig,ax = subplots(length(hists), 1, figsize=figsize)
        length(hists) == 1 && (ax = [ax])
        for (i, h) = enumerate(hists)
            bwidth = diff(h.edges[1])[1]
            bx = h.edges[1][1:end-1] .+ bwidth/2
            ax[i].bar(bx, h.weights, width=bwidth, edgecolor="black")
        end
        nicenup(fig, fsize=fsize)
    end
    hists, CI
end

function makenuisancehists(optn::OptionsNuisance, qp1, qp2; burninfrac = 0.5, nbins = 50,
    temperaturenum = -1, pdfnormalize = false)
    @assert temperaturenum != -1 "Please explicitly specify the temperature index"
    nuisanceatT = assemblenuisancesatT(optn, burninfrac = burninfrac,
                                    temperaturenum = temperaturenum)
    nnu = length(optn.idxnotzero)
    CI = zeros(Float64, nnu, 3)
    h = Vector{Histogram}(undef, nnu)
    for (i, islice) in enumerate(optn.idxnotzero)
        numin, numax = extrema(optn.bounds[islice,:])
        edges = LinRange(numin, numax, nbins+1)
        h[i] = fit(Histogram, nuisanceatT[:,islice], edges)
        CI[i,:] = quantile(nuisanceatT[:,islice], [qp1, 0.5, qp2])
    end
    h, CI
end

function firstderiv(x)
    fd = similar(x)
    fd[1]        =  x[2]     - x[1]
    fd[end]      =  x[end]   - x[end-1]
    fd[2:end-1] .= (x[3:end] - x[1:end-2])/2
    abs.(fd)
    # fd
end

function secondderiv(x)
    sd = similar(x)
    sd[1]        = 2x[1] - 5x[2] + 4x[3] - x[4]
    sd[2:end-1] .= (x[3:end] - x[2:end-1] + x[1:end-2])
    sd[end]      = 2x[end] - 5x[end-1] + 4x[end-2] - x[end-3]
    abs.(sd)
end

function gethimage(F::Operator, M::AbstractArray, opt::Options;
                burninfrac = 0.5,
                nbins = 50,
                rhomin=Inf,
                rhomax=-Inf,
                qp1=0.05,
                qp2=0.95,
                islscale = false,
                pdfnormalize=false,
                istothepow = false,
                temperaturenum=1)
    T = x->x
    if (rhomin == Inf) && (rhomax == -Inf)
        if !stretchexists(F) # if no stretch
            for (i,mm) in enumerate(M)
                rhomin_mm = minimum(mm)
                rhomax_mm = maximum(mm)
                rhomin_mm < rhomin && (rhomin = rhomin_mm)
                rhomax_mm > rhomax && (rhomax = rhomax_mm)
            end
            if (typeof(opt) == OptionsStat && opt.needλ²fromlog) && islscale
                T = x->0.5log10(x)
            end
        else # there is a stretch
            rhomin = minimum(F.low)
            rhomax = maximum(F.low + F.Δ)
        end        
    else
        @assert rhomin < rhomax
    end
    if istothepow && !islscale
        T = x->10. ^x
    end
    rhomin, rhomax = map(x->T(x), (rhomin, rhomax))    
    edges = LinRange(rhomin, rhomax, nbins+1)
    himage = zeros(Float64, length(M[1]), nbins)
    CI = zeros(Float64, length(M[1]), 3)
    meanimage = zeros(length(M[1]))
    for ilayer=1:length(M[1])
        if !stretchexists(F) # if no stretch
            mthislayer = [T(m[ilayer]) for m in M]
            himage[ilayer,:] = fit(Histogram, mthislayer, edges).weights
            CI[ilayer,:] = [quantile(mthislayer,(qp1, 0.5, qp2))...]
            meanimage[ilayer] = mean(vec(mthislayer))
        else # there is an affine stretch with depth
            mthislayer = [m[ilayer] for m in M]
            expandedthislayer = T.(F.low[ilayer] .+ mthislayer*F.Δ[ilayer])
            himage[ilayer,:] = fit(Histogram, expandedthislayer, edges).weights
            CI[ilayer,:] = [quantile(expandedthislayer,(qp1, 0.5, qp2))...]
            meanimage[ilayer] = mean(vec(expandedthislayer))
        end        
        himage[ilayer,:] = himage[ilayer,:]/sum(himage[ilayer,:])/(diff(edges)[1])
        pdfnormalize && (himage[ilayer,:] = himage[ilayer,:]/maximum(himage[ilayer,:]))
    end
    if !stretchexists(F)
        Tm = [T.(m) for m in M]
    else
        Tm = [T.(F.low + vec(m).*F.Δ) for m in M]
    end    
    Mslope = mean(firstderiv.(Tm))
    sdevslope = std(firstderiv.(Tm))
    himage, edges, CI, meanimage, Mslope, sdevslope
end

function stretchexists(F::Operator)
    in(:stretch, fieldnames(typeof(F))) && F.stretch
end

# plotting codes for 2D sections in AEM

function splitsoundingsbyline(soundings::Array{S, 1}) where S<:Sounding
    alllines = [s.linenum for s in soundings]
    linenum  = unique(alllines)
    nlines   = length(linenum)
    linestartidx = zeros(Int, nlines)
    for i = 1:nlines
        linestartidx[i] = findfirst(alllines .== linenum[i])
    end
    linestartidx
end 

function makegrid(vals::AbstractArray, soundings::Array{S, 1};
    dr=10, zall=[NaN], dz=-1) where S<:Sounding
    @assert all(.!isnan.(zall)) 
    @assert dz>0
    X = [s.X for s in soundings]
    Y = [s.Y for s in soundings]
    x0, y0 = X[1], Y[1]
    R  = sqrt.((X .- x0).^2 + (Y .- y0).^2)
    rr, zz = [r for z in zall, r in R], [z for z in zall, r in R]
    topo = [s.Z for s in soundings]
    zz = topo' .- zz # mAHD
    kdtree = KDTree([rr[:]'; zz[:]'])
    gridr = range(R[1], R[end], step=dr)
    gridz = reverse(range(extrema(zz)..., step=dz))
    rr, zz = [r for z in gridz, r in gridr], [z for z in gridz, r in gridr]
    idxs, = nn(kdtree, [rr[:]'; zz[:]'])
    img = zeros(size(rr))
    for i = 1:length(img)
        img[i] = vals[idxs[i]]
    end
    topofine = gridpoints(R, gridr, topo)
    img[zz .>topofine'] .= NaN
    img, gridr, gridz, topofine, R
end

function makegrid(vals::AbstractArray, soundings::Array{S, 1}, XYZrho;
    dr=10, zall=[NaN], dz=-1) where S<:Sounding
    @assert all(.!isnan.(zall)) 
    @assert dz>0
    X = [s.X for s in soundings]
    Y = [s.Y for s in soundings]
    x0, y0 = X[1], Y[1]
    R  = sqrt.((X .- x0).^2 + (Y .- y0).^2)
    rr, zz = [r for z in zall, r in R], [z for z in zall, r in R]
    topo = [s.Z for s in soundings]
    zz = topo' .- zz # mAHD
    kdtree = KDTree([rr[:]'; zz[:]'])
    gridr = range(R[1], R[end], step=dr)
    gridz = reverse(range(extrema(zz)..., step=dz))
    rr, zz = [r for z in gridz, r in gridr], [z for z in gridz, r in gridr]
    idxs, = nn(kdtree, [rr[:]'; zz[:]'])
    # comparison with ZYZrho
    x0, y0 = XYZrho[1,1], XYZrho[1,2]
    rrcompare = sqrt.((XYZrho[:,1] .- x0).^2 + (XYZrho[:,2] .- y0).^2)
    kdtreecompare = KDTree([rrcompare[:]'; XYZrho[:,3]'])
    idxscompare, = nn(kdtreecompare, [rr[:]'; zz[:]'])
    img, imgcompare = zeros(size(rr)), zeros(size(rr))
    for i = 1:length(img)
        img[i] = vals[idxs[i]]
        imgcompare[i] = XYZrho[idxscompare[i], 4]
    end
    topofine = gridpoints(R, gridr, topo)
    img[zz .>topofine'] .= NaN
    imgcompare[zz .> topofine'] .= NaN
    img, gridr, gridz, topofine, R, imgcompare
end

function gridpoints(R, gridr, points::Array{Float64, 1})
    @assert length(R) == length(points)
    kdtree = KDTree(R[:]')
    idxs, = nn(kdtree, gridr[:]')
    pointsgrid = [points[idxs[i]] for i = 1:length(gridr)]
end

function gridpoints(R, gridr, points::Array{Float64, 2})
    nsound = length(R)
    nprops = size(points, 1)
    @assert nsound == size(points, 2)
    pointsgrid = zeros(length(gridr), nprops)
    for i = 1:nprops
        pointsgrid[:,i] = gridpoints(R, gridr, vec(points[i,:]))
    end
    pointsgrid    
end    

function whichislast(soundings::AbstractArray)
    X = [s.X for s in soundings]
    Y = [s.Y for s in soundings]
    Eislast, Nislast = true, true
    X[1]>X[end] && (Eislast = false)
    Y[1]>Y[end] && (Nislast = false)
    Eislast, Nislast
end

function makesummarygrid(soundings, pl, pm, ph, ρmean, vdmean, vddev, zall, dz; dr=10)
    # first flip ρ to σ and so pl and ph are interchanged
    phgrid, gridx, gridz, topofine, R = makegrid(-pl, soundings, zall=zall, dz=dz, dr=dr)
    plgrid,                           = makegrid(-ph, soundings, zall=zall, dz=dz, dr=dr)
    pmgrid,                           = makegrid(-pm, soundings, zall=zall, dz=dz, dr=dr)
    σmeangrid,                        = makegrid(-ρmean, soundings, zall=zall, dz=dz, dr=dr)
    ∇zmeangrid,                       = makegrid(vdmean, soundings, zall=zall, dz=dz, dr=dr)
    ∇zsdgrid,                         = makegrid(vddev, soundings, zall=zall, dz=dz, dr=dr)
    Z = [s.Z for s in soundings]
    phgrid, plgrid, pmgrid, σmeangrid, ∇zmeangrid, ∇zsdgrid, gridx, gridz, topofine, R, Z
end

function makexyzrho(soundings, zall)
    linename = "_line_$(soundings[1].linenum)_summary.txt"
    fnames = ["rho_low", "rho_mid", "rho_hi", "rho_avg"].*linename
    dlow, dmid, dhigh, davg = map(x->readdlm(x), fnames)
    ndepths = length(zall)
    for (i, d) in enumerate([dlow, dmid, dhigh, davg])
        @assert ndims(d) == 2
        @assert size(d, 2) == length(soundings)
        @assert size(d, 1) == ndepths
        xyzrho = makearray(soundings, d, zall)
        writedlm(fnames[i][1:end-4]*"_xyzrho.txt", xyzrho)
    end
end

function makearray(soundings, d, zall)
    outarray = zeros(length(d), 4)
    ndepths = length(zall)
        for (i, s) in enumerate(soundings)
            z = s.Z .- zall # convert depths to mAHD
            x, y = s.X, s.Y
            idx = (i-1)*ndepths+1:i*ndepths
            outarray[idx,:] = [[x y].*ones(ndepths) z d[:,i]]
        end
    outarray
end

function plotNEWSlabels(Eislast, Nislast, gridx, gridz, axarray)
    for s in axarray
        Eislast ? s.text(gridx[1], minimum(s.get_ylim()), "W", backgroundcolor="w") : s.text(gridx[1], minimum(s.get_ylim()), "E", backgroundcolor="w")
        Nislast ? s.text(gridx[end], minimum(s.get_ylim()), "N", backgroundcolor="w") : s.text(gridx[end], minimum(s.get_ylim()), "S", backgroundcolor="w")
    end
end

function plotprofile(ax, idxs, Z, R)
    for idx in idxs
        ax.plot(R[idx]*[1,1], [ax.get_ylim()[1], Z[idx]], "-w")
        ax.plot(R[idx]*[1,1], [ax.get_ylim()[1], Z[idx]], "--k")
    end
end

function selectwithin1Dinterval(M::AbstractVector, z, 
                                zbounds, vbounds,
                                cond = :mean)
    @assert(any(cond .== [:all, :any, :median, :mean]))
    @assert size(vbounds, 1) == size(zbounds, 1)
    @assert length(z) == length(M[1])
    @assert all(vbounds[:,1] .< vbounds[:,2])
    @assert all(zbounds[:,1] .< zbounds[:,2])
    nconditions = size(vbounds, 1)
    idx = zeros(Bool, length(M), nconditions)
    for (i, m) in enumerate(M)
       for j in 1:nconditions
            idxdepth = zbounds[j,1] .< z .< zbounds[j,2]
            if any(cond .== [:median, :mean])
                if vbounds[j,1] < eval(cond)(m[idxdepth]) < vbounds[j,2]
                    idx[i,j] = true
                end
            else     
                if eval(cond)(vbounds[j,1] .< (m[idxdepth]) .< vbounds[j,2])
                    idx[i,j] = true
                end    
            end
       end
    end
    vec(reduce(&, idx, dims=2))
end

function block1Dvalues(M::AbstractVector, z, zbounds, cond = :median)
    @assert(any(cond .== [:median, :mean]))
    @assert length(z) == length(M[1])
    @assert all(zbounds[:,1] .< zbounds[:,2])
    nconditions = size(zbounds, 1)
    Mblock = zeros(length(M), nconditions)
    for (i, m) in enumerate(M)
        for j in 1:nconditions
            idxdepth = zbounds[j,1] .< z .<= zbounds[j,2]
            Mblock[i,j] = eval(cond)(m[idxdepth])
        end
     end
    Mblock
end    

function correlationplot(M::Array{Float64, 2}; figsize=(5,5), nbins=25)
    f, ax = plt.subplots(size(M, 2), size(M, 2), figsize=figsize, sharex=true, sharey=true)
    nrows = size(M, 2)
    for j = 1:nrows
        for i = j:nrows
            h = normalize(fit(Histogram, (M[:,i], M[:,j]), nbins=nbins))
            ax[nrows*(j-1)+i].pcolormesh(h.edges[2], h.edges[1], h.weights)
        end
    end 
    ax[1].set_aspect(1)   
end    

function plot1Dpatches(ax, zlist, xlist; alpha=0.25, fc="red", ec="blue", lw=2)
    @assert length(zlist) == length(xlist)
    for i in 1:size(zlist, 1)
        x0, z0 = xlist[i,1], zlist[i,1]
        delx = xlist[i,2] - x0
        delz = zlist[i,2] - z0
        ax.add_patch(matplotlib.patches.Rectangle((x0,z0),delx,delz, alpha=alpha, fc=fc, ec=fc))
        ax.add_patch(matplotlib.patches.Rectangle((x0,z0),delx,delz, fill=false, ec=ec, lw=lw, linestyle="--"))
    end    
end

#filtering the projection added part from oasis 
function filter_proj!(df)
    filter!(row -> !contains(row.combined, "RT=PROJ"), df)
end

#filtering the transformation added part from oasis 
function filter_trans!(df)
    filter!(row -> !contains(row.combined,"RT=TRNS"), df)
end 

#function to read the *dfn file and extract the column number and column names as a *.txt file 
function dfn2hdr(dfnfile::String)
    df = CSV.read(dfnfile, DataFrame; header=false);
    df = coalesce.(df,"")
    
    colname = []
    for col in 1:length(names(df))
        push!(colname,names(df)[col])
    end 
    t = df[!, colname[1]]
    for i in colname[2:end]         
		t = t .* df[!, i]                  
	end 
    df[!,"combined"] = t
    filter_proj!(df)
    filter_trans!(df)
    row_number = nrow(df)
    insertcols!(df, 1, :inc => 0)
    insertcols!(df, 1, :first => 0)
    insertcols!(df,1, :second => 0)

    cumulative_columns = 0  #this will set up a cumulative variable 
    
    for i in 2:size(df,1)   #this will do operations in row level 
        regex_literal1 = r"NAME=" #few parameters have their names after this keyword 
        regex_literal2 = r"RT="  #four parameters have their names after this keyword 
        if !contains(string(df[i,:].combined),regex_literal1)
            continue
        end

        if contains(df[i,:].combined,"RT:A4F")
            df[i,:].combined = replace(df[i,:].combined, ":A4F" => "AA")
        end 

        inc_regex = r":([0-9]+)F"
        inc_match = match(inc_regex, df[i,:].combined)  #matching the keyword with the string 
        if isnothing(inc_match)
            inc = 0
            idx2 =  idx2 = split(string(df[i,:].combined), regex_literal1)
            idx2_r = first.(split.(idx2[2], ":"))
            
        else
            inc = parse(Int64, inc_match.captures[1]) #the outcome is int64
            idx2 =  idx2 = split(string(df[i,:].combined), regex_literal2)
            idx2_r = first.(split.(idx2[2], ":"))
           
        end

        df[i,:].inc = inc
        df[i,:].first = cumulative_columns + 1
        df[i,:].second = df[i,:].first + inc 
        cumulative_columns += inc + 1
        
        if occursin(";END DEFN", idx2_r)  #type is string 
            idx2_r = first.(split.(idx2_r,";"))
        else
            idx2_r = first.(split.(idx2_r,"END"))
        end

        #start writing to a file here 
        io = open("header.txt","a")
        if (df[i,:].first != df[i,:].second)
            writedlm(io,[df[i,:].first df[i,:].second idx2_r]) 
        else
            writedlm(io, [df[i,:].first idx2_r])
        end
        close(io)
    end
end

end # module CommonToAll