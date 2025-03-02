module CommonToAll
using PyPlot, StatsBase, Statistics, Distances, LinearAlgebra,
      DelimitedFiles, ..AbstractOperator, NearestNeighbors, Printf, ReadVTK,
      KernelDensitySJ, KernelDensity, Interpolations, CSV, WriteVTK, Distributed, Glob

import ..Options, ..OptionsStat, ..OptionsNonstat, ..OptionsNuisance,
       ..history, ..GP.κ, ..calcfstar!, ..AbstractOperator.Sounding, 
       ..AbstractOperator.getsmoothline,
       ..DEBUGLEVEL_TDGP

export trimxft, assembleTat1, gettargtemps, checkns, getchi2forall, nicenup, plotconv,
        plot_posterior, make1Dhist, make1Dhists, setupz, zcontinue, makezρ, plotdepthtransforms,
        unwrap, getn, geomprogdepth, assemblemodelsatT, getstats, gethimage,
        assemblenuisancesatT, makenuisancehists, stretchexists, stepmodel,
        makegrid, whichislast, makesummarygrid, makearray, plotNEWSlabels, 
        plotprofile, gridpoints, splitsoundingsbyline, getsoundingsperline, docontinue, linestartend,
        compatidxwarn, dfn2hdr, getgdfprefix, readlargetextmatrix, pairinteractionplot, flipline, 
        summaryconductivity, plotsummarygrids1, getVE, writevtkfromsounding, trapezoidrule,
        readcols, colstovtk, findclosestidxincolfile, zcentertoboundary, zboundarytocenter, zboundarytocenter_inexact,
        writeijkfromsounding, nanmean, infmean, nanstd, infstd, infnanmean, infnanstd, calchighalt,
        kde_sj, plotmanygrids, readwell, getlidarheight, plotblockedwellonimages, getdeterministicoutputs, writeasegdfnfromonerow,
        writeasegdat, getprobabilisticoutputs, readfzipped, readxyzrhoϕ, writevtkxmlforcurtain, getRandgridr, getallxyinr, getXYlast,
        getprobabilisticlinesfromdirectory, readxyzrhoϕnu, plotgausshist, findMAE, getclosestidx, getpostatXY, getpostatwell, thicktodepth

# Kernel Density stuff
abstract type KDEtype end
struct SJ <: KDEtype end
struct LSCV <: KDEtype end

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
    opt.fstar_filename = "models_"*opt.fdataname*isns*".bin"
    opt.x_ftrain_filename = "points_"*opt.fdataname*isns*".bin"
    opt.costs_filename = "misfits_"*opt.fdataname*isns*".bin"
    chain_idx = nothing
    if isfile(opt.fstar_filename)
        chain_idx = 1
    end
    for ichain in 1:size(Tacrosschains, 2)
        if isnothing(chain_idx)
            opt.fstar_filename = "models_"*opt.fdataname*isns*"_$ichain.bin"
            opt.x_ftrain_filename = "points_"*opt.fdataname*isns*"_$ichain.bin"
            opt.costs_filename = "misfits_"*opt.fdataname*isns*"_$ichain.bin"
        else
            chain_idx = ichain
        end
        if stat == :fstar
            at1idx = findall(Tacrosschains[:,ichain].==ttarget) .>= start
        else
            at1idx = Tacrosschains[:,ichain].==ttarget
            at1idx[1:start-1] .= false
        end
        ninchain = sum(at1idx)
        (DEBUGLEVEL_TDGP > 0) && @info("chain $ichain has $ninchain models")
        ninchain == 0 && continue
        mat1[imodel+1:imodel+ninchain] .= history(opt, stat=stat, chain_idx=chain_idx)[at1idx]
        imodel += ninchain
    end
    iters = history(opt, stat=:iter, chain_idx=chain_idx)
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
    firsti == 0 && (firsti = 1)
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
        (DEBUGLEVEL_TDGP > 0) && @info("chain $ichain has $ninchain models")
        ninchain == 0 && continue
        vals_filename = "values_nuisance_"*fdataname*".bin"
        if isfile(vals_filename)
            nraw = readdlm(vals_filename, ' ', String)
            cids = parse.(Int, nraw[:,1])
            ndat = parse.(Float64, nraw[cids .== ichain, 3:end])
        else
            vals_filename = "values_nuisance_"*fdataname*"$ichain.bin"
            ndat = readdlm(vals_filename, ' ', Float64)[:,2:end]
        end
        mvals[imodel+1:imodel+ninchain,:] .= ndat[at1idx,:]
        imodel += ninchain
    end
    mvals
end

function getnchains(costs_filename)
    if isfile(costs_filename*".bin")
        data = readdlm(costs_filename * ".bin", String)
        chids = parse.(Int, data[:,1])
        return maximum(chids)
    end

    c = 0
    r = Regex(costs_filename)
    for fname in readdir(pwd())
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
    multichainfile = nothing
    if isfile(costs_filename*".bin")
        opt.costs_filename = costs_filename*".bin"
        multichainfile = 1 # set if all chains are in one file
    else
        opt.costs_filename    = costs_filename*"_1.bin"
    end
    iters          = history(opt, stat=:iter, chain_idx=multichainfile)
    niters         = length(iters)
    @info "McMC has run for $(iters[end]) iterations"
    # then create arrays of unsorted by temperature T
    Tacrosschains  = zeros(Float64, niters, nchains)
    # get the values into the arrays
    for ichain in 1:nchains
        if isnothing(multichainfile)
            opt.costs_filename = costs_filename*"_$ichain.bin"
            Tacrosschains[:,ichain] = history(opt, stat=:T)
        else
            Tacrosschains[:,ichain] = history(opt, stat=:T, chain_idx=ichain)
        end
    end
    Tacrosschains
end

function gettargtemps(optn_in::OptionsNuisance)
    optn = deepcopy(optn_in)
    costs_filename = "misfits_"*optn.fdataname*"nuisance"
    nchains = getnchains(costs_filename)
    @info "Number of chains is $nchains"
    multichainfile = nothing
    if isfile(costs_filename*".bin")
        optn.costs_filename = costs_filename*".bin"
        multichainfile = 1 # set if all chains are in one file
    else
        optn.costs_filename    = costs_filename*"_1.bin"
    end
    iters = history(optn, stat=:iter, chain_idx=multichainfile)
    niters = length(iters)
    @info "MCMC has run for $(iters[end]) iterations"
    Tacrosschains = zeros(Float64, niters, nchains)
    for ichain in 1:nchains
        if isnothing(multichainfile)
            optn.costs_filename = costs_filename*"_$ichain.bin"
            Tacrosschains[:,ichain] = history(optn, stat=:T)
        else
            Tacrosschains[:,ichain] = history(optn, stat=:T, chain_idx=ichain)
        end
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
    multichainfile = nothing
    if isfile(costs_filename*".bin")
        opt.costs_filename = costs_filename*".bin"
        multichainfile = 1
    else
        opt.costs_filename    = costs_filename*"_1.bin"
    end
    iters          = history(opt, stat=:iter, chain_idx=multichainfile)
    statnames = [:acceptanceRateBirth, :acceptanceRateDeath,
                 :acceptanceRatePosition, :acceptanceRateProperty, :acceptanceRateDC]
    f,ax = plt.subplots(length(statnames), 1,
                        sharex=true, sharey=true, figsize=figsize)
    maxar = 0
    for (ichain, chain) in enumerate(chains)
        if isnothing(multichainfile)
            opt.costs_filename = costs_filename*"_$chain.bin"
            chain_idx = nothing
        else
            chain_idx = ichain
        end
        for (istat, statname) in enumerate(statnames)
            ar = history(opt, stat=statname, chain_idx=chain_idx)
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
    multichainfile = nothing
    if isfile(costs_filename*".bin")
        opt.costs_filename = costs_filename*".bin"
        multichainfile = 1
    else
        opt.costs_filename = costs_filename*"_1.bin"
    end
    iters          = history(opt, stat=:iter, chain_idx=multichainfile)
    statname = :acceptanceRate
    f,ax = plt.subplots(length(optin.idxnotzero), 1, sharex=true, sharey=true, figsize=figsize)
    maxar = 0
    for (ichain, chain) in enumerate(chains)
        if isnothing(multichainfile)
            opt.costs_filename = costs_filename*"_$chain.bin"
            chain_idx = nothing
        else
            chain_idx = ichain
        end
        ar = history(opt, stat=statname, chain_idx=chain_idx)[:,optin.idxnotzero]
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
    (DEBUGLEVEL_TDGP > 0) && @info("ns is $isns")
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
    multichainfile = false
    if isfile(optn.costs_filename)
        iters = history(optn, stat=:iter)
    else
        multichainfile = true
        optn.costs_filename = costs_filename*".bin"
        iters = history(optn, stat=:iter, chain_idx=1)
    end
    niters = length(iters)
    Tacrosschains = zeros(Float64, niters, nchains)
    χ2acrosschains = zeros(Float64, niters, nchains)
    chain_idx = nothing
    for chain in 1:nchains
        if multichainfile
            chain_idx = ichain
        else
            optn.costs_filename = costs_filename*"_$(chain).bin"
        end
        Tacrosschains[:,chain] = history(optn, stat=:T, chain_idx=chain_idx)
        χ2acrosschains[:,chain] = history(optn, stat=:misfit, chain_idx=chain_idx)
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
                        alpha           = 0.5,
                        nxticks         = 0,
                        gridon          = false,
                        omittemp        = false,
                        hidetitle       = true,
                      )
    # now look at any chain to get how many iterations
    isns = checkns(opt_in)
    opt = deepcopy(opt_in)
    costs_filename = "misfits_"*opt.fdataname*isns
    if nchains == 1 # then actually find out how many chains there are saved
        nchains = getnchains(costs_filename)
    end
    opt.costs_filename    = costs_filename*"_1.bin"
    multichainfile = false
    if isfile(opt.costs_filename)
        iters          = history(opt, stat=:iter)
    else
        opt.costs_filename = costs_filename*".bin"
        iters = history(opt, stat=:iter, chain_idx=1)
        multichainfile=true
    end
    niters         = length(iters)
    # then create arrays of unsorted by temperature T, k, and chi2
    Tacrosschains  = zeros(Float64, niters, nchains)
    kacrosschains  = zeros(Int, niters, nchains)
    X2by2inchains  = zeros(Float64, niters, nchains)
    # get the values into the arrays
    chain_idx = nothing
    for ichain in 1:nchains
        if multichainfile
            chain_idx = ichain
        else
            opt.costs_filename = costs_filename*"_$ichain.bin"
        end
        Tacrosschains[:,ichain] = history(opt, stat=:T, chain_idx=chain_idx)
        kacrosschains[:,ichain] = history(opt, stat=:nodes, chain_idx=chain_idx)
        X2by2inchains[:,ichain] = history(opt, stat=:U, chain_idx=chain_idx)
    end

    if !hidetitle # then we are usually not interested in the temperature sorting of chains    
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
    end

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
    !hidetitle && ax[1].set_title(isns*" sorted by temperature")
    ax[1].plot(iters, kacrosschains[:,1:nchainsatone], "k", alpha=alpha)
    gridon && ax[1].grid()
    chi2low, chi2high = 0.01*median(X2by2inchains[:,1]), 5*median(X2by2inchains[:,end])
    ax[2].plot(iters, X2by2inchains, alpha=alpha)
    ax[2].plot(iters, X2by2inchains[:,1:nchainsatone], "k", alpha=alpha)
    ax[2].set_ylabel("-Log L")
	ax[2].set_ylim(chi2low, chi2high)
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

function plotconv(optin::Options; burninfrac= 0.5, first=0.1, last=0.5, till=1.0, figsize=(6,2), fontsize=10, nbins=50, doall=false)
    opt = deepcopy(optin)
    iter, k, chi2by2 = map(x->assembleTat1(opt, x; temperaturenum=1, burninfrac), (:iter, :nodes, :U))
    idx = sortperm(iter)
    till = round(Int, idx[end]*till)
    iter, k, chi2by2 = map(x->x[idx][1:till], (iter, k, chi2by2))
    f, ax = plt.subplots(1, 2; figsize)
    # first group, second group
    getconv(ax, iter, chi2by2, k, nbins, opt.nmin, opt.nmax, first, last)
    ax[1].set_xlabel("-ve log likelihood")
    ax[1].set_ylabel("pdf")
    ax[2].set_xlabel("# nuclei")
    ax[2].set_ylabel("pdf")
    nicenup(f, fsize=fontsize)
    nothing
end 

function getconv(ax, iter, chi2by2, k, nbins, nmin, nmax, first, last)
    first, last = map(x->round(Int, length(iter)*x), (first, last))
    plotconv(ax, chi2by2[1:first], k[1:first], nbins, nmin, nmax)
    plotconv(ax, chi2by2[last:end], k[last:end], nbins, nmin, nmax)
end    

function plotconv(ax, chi2by2, k, nbins, nmin, nmax)
    kdefunc_nll = kde_(SJ(), chi2by2)
    nllrange = range(minimum(chi2by2), maximum(chi2by2), nbins)
    krange = nmin-0.5:1:nmax+0.5
    ax[1].plot(nllrange, kdefunc_nll(nllrange))
    ax[2].hist(k, krange, density=true, histtype="step")
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

function plotdepthtransforms(zall, znall, zboundaries; fontsize=12)
    thickness = diff(zboundaries)
    f, ax = plt.subplots(2, 2, sharex="col", sharey="row", figsize=(8,5))
    ax[1].stem(zboundaries[1:end-1], zboundaries[1:end-1], markerfmt="")
    ax[1].stem(zall, zall, "k--", markerfmt=" ")
    ax[1].set_ylabel("depth m")
    ax[2].stem(zall,thickness, "k--", markerfmt=" ")
    ax[2].set_ylabel("thickness m")
    ax[2].set_xlabel("depth m")
    ax[2].yaxis.grid(which="major")
    ax[3].plot(znall, zall)
    ax[3].grid()
    ax[4].stem(znall,thickness, "k--", markerfmt=" ")
    ax[4].set_xlabel("depth index")
    ax[4].yaxis.grid(which="major")
    plt.suptitle("Forward model discretization"; fontsize)
    nicenup(f, fsize=fontsize)
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

function nicenup(g::PyPlot.Figure;fsize=12, h_pad=nothing, w_pad = nothing, increasefraction=1.1, minsize=true)
    for ax in g.axes
        if !isempty(ax.get_yticklabels())
            fs = ax.get_yticklabels()[1].get_fontsize()
            ns = getnewfontsize(fs, increasefraction, fsize; minsize)
            ax.tick_params("both", labelsize=ns)
        elseif !isempty(ax.get_xticklabels())
            fs = ax.get_xticklabels()[1].get_fontsize()
            ns = getnewfontsize(fs, increasefraction, fsize; minsize)
            ax.tick_params("both", labelsize=ns)
        end
        xlh, ylh, tlh = ax.xaxis.label, ax.yaxis.label, ax.title
        for h in (xlh, ylh, tlh)
            fs = h.get_fontsize()
            ns = getnewfontsize(fs, increasefraction, fsize; minsize)
            h.set_fontsize(ns)
        end    
        if any(keys(ax) .== :zaxis)
            fs = ax.zaxis.label.get_fontsize()
            ns = getnewfontsize(fs, increasefraction, fsize; minsize)
            ax.zaxis.label.set_fontsize(ns)
        end    
        if typeof(ax.get_legend_handles_labels()[1]) != Array{Any,1}
            ax.legend(loc="best", fontsize=fsize)
            fs = ax.get_legend().get_texts()[1].get_fontsize()
            ns = getnewfontsize(fs, increasefraction, fsize; minsize)
            ax.legend(loc="best", fontsize=ns)
        end
    end
    if isnothing(h_pad)
        g.tight_layout()
    else
        g.tight_layout(;h_pad)
    end
    if isnothing(w_pad)
        g.tight_layout()
    else
        g.tight_layout(;w_pad)
    end 
    nothing
end

function getnewfontsize(fs, increasefraction, fsize; minsize=true)
    if minsize
        ns = fs*increasefraction
        ns = ns > fsize ? ns : fsize
    else #exactsize
        ns = fsize
    end
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
                        figsize=(10,5),
                        pdfnormalize=false,
                        fsize=14,
                        istothepow=false,
                        usekde = false,
                        kdetype = SJ(),
                        alpha=0.25,
                        lw = 1)
    @assert 0<vmaxpc<=1
    M = assembleTat1(optns, :fstar, burninfrac=burninfrac, temperaturenum=temperaturenum)
    himage_ns, edges_ns, CI_ns, = gethimage(F, M, optns; temperaturenum=temperaturenum,
                            nbins=nbins, qp1=qp1, qp2=qp2, istothepow=istothepow, usekde, kdetype,
                            islscale=false, pdfnormalize=pdfnormalize)
    M = assembleTat1(opts, :fstar, burninfrac=burninfrac, temperaturenum=temperaturenum)
    himage, edges, CI,  = gethimage(F, M, opts; temperaturenum=temperaturenum,
                            nbins=nbins, qp1=qp1, qp2=qp2, istothepow=false, usekde, kdetype,
                            islscale=true, pdfnormalize=pdfnormalize)
    f,ax = plt.subplots(1, 2, sharey=true, figsize=figsize)
    xall = opts.xall
    xmesh = [zcentertoboundary(xall); xall[end]]
    vmin, vmax = extrema(himage_ns)
    vmax = vmin+vmaxpc*(vmax-vmin)
    im1 = ax[1].pcolormesh(edges_ns[:], xmesh, himage_ns, cmap=cmappdf, vmax=vmax)
    cb1 = colorbar(im1, ax=ax[1])
    cb1.ax.set_xlabel("pdf \nns")
    propmin, propmax = getbounds(CI_ns, optns.fbounds)
    ax[1].set_xlim(propmin, propmax)
    
    vmin, vmax = extrema(himage)
    vmax = vmin+vmaxpc*(vmax-vmin)
    im2 = ax[2].pcolormesh(edges[:], xmesh, himage, cmap=cmappdf, vmax=vmax)
    ax[2].set_ylim(extrema(xall)...)
    propmin, propmax = getbounds(CI, opts.fbounds)
    ax[2].set_xlim(propmin, propmax)
    ax[2].invert_yaxis()
    cb2 = colorbar(im2, ax=ax[2])
    cb2.ax.set_xlabel("pdf \nstationary")
    ax[1].plot(CI_ns, xall[:], linewidth=lw, color="w"; alpha)
    ax[2].plot(CI, xall[:], linewidth=lw, color="w"; alpha)
    ax[1].set_xlabel(L"\log_{10} \rho")
    ax[1].set_ylabel("depth (m)")
    ax[2].set_xlabel(L"\log_{10} λ")
    nicenup(f, fsize=fsize)
end

function plot_posterior(F::Operator1D,
                    opt::OptionsStat;
                    useoptfbounds = true,
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
                    pdfclim = nothing,
                    showslope = false,
                    kdetype = SJ(),
                    usekde = false,
                    doplot = true)
    @assert 0<vmaxpc<=1
    
    if useoptfbounds
        if stretchexists(F)
            rhomin, rhomax = minimum(F.low), maximum(F.low + F.Δ)
        else
            rhomin, rhomax = extrema(opt.fbounds)
        end
    else
        rhomin, rhomax = Inf, -Inf
    end    
    M = assembleTat1(opt, :fstar, burninfrac=burninfrac, temperaturenum=temperaturenum)
    himage, edges, CI, meanimage, meandiffimage, sdslope, = gethimage(F, M, opt; temperaturenum,
                nbins, qp1, qp2, istothepow, rhomin, rhomax, usekde, kdetype,
                islscale=false, pdfnormalize=pdfnormalize)

    if doplot
        if showslope
            f, ax = plt.subplots(1,2, sharey=true, figsize=figsize)
        else
            f, ax = plt.subplots(1,1, sharey=true, figsize=figsize, squeeze=false)
        end    
        xall = opt.xall
        xmesh = [zcentertoboundary(xall); xall[end]]
        vmin, vmax = extrema(himage)
        vmax = vmin+vmaxpc*(vmax-vmin)
        im1 = ax[1].pcolormesh(edges[:], xmesh, himage, cmap=cmappdf, vmax=vmax)
        ax[1].set_ylim(extrema(xall)...)
        ax[1].invert_yaxis()
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
        cb1 = colorbar(im1, ax=ax[1])
        cb1.ax.set_title("pdf")
        if showslope
            ax[2].plot(meandiffimage[:], xall[:], linewidth=2, color="k", linestyle="-")
            zeroside = meandiffimage[:]-sdslope[:]
            zeroside[zeroside .< 0] .= 0
            ax[2].fill_betweenx(xall[:],zeroside,meandiffimage[:]+sdslope[:],alpha=.25)
            ax[2].set_xlabel("mean slope")
        end
        !isnothing(pdfclim) && ax[1].collections[1].set_clim(pdfclim)
        nicenup(f, fsize=fsize)
    end
    CI[:,1], CI[:,2], CI[:,3], meanimage, meandiffimage, sdslope
end

function getpostatXY(F::Operator, soundings::Vector{T}, opt_in::OptionsStat, Xwanted, Ywanted; 
        burninfrac = 0.5,
        nbins = 50, 
        qp1 = nothing, 
        qp2 = nothing,
        usekde = true,
        pdfnormalize = false,
        ) where T<:Sounding
    @assert !isnothing(qp1) | !isnothing(qp2)
    opt = deepcopy(opt_in)
    M = getclosestensemble(soundings, opt, Xwanted, Ywanted; burninfrac)
    rhomin, rhomax = extrema(opt.fbounds)
    himage, edges, CI, _ = gethimage(F, M, opt; temperaturenum=1,
                nbins, qp1, qp2, rhomin, rhomax, usekde,
                islscale=false, pdfnormalize)
    himage, edges, CI
end

function getpostatwell(F::Operator, soundings::Vector{T}, opt_in::OptionsStat, Xwell, Ywell, zc_rho; 
    burninfrac = 0.5,
    nbins = 50, 
    qp1 = nothing, 
    qp2 = nothing,
    usekde = true,
    pdfnormalize = false,
    ) where T<:Sounding
    zcentre = opt_in.xall
    Mwell = blocktomodel(zcentre, zc_rho)
    himage, edges, CI = getpostatXY(F, soundings, opt_in, Xwell, Ywell; burninfrac,
                                nbins, qp1, qp2, usekde, pdfnormalize)
    himage, edges, CI, Mwell
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
                        figsize=(5,8),
                        pdfnormalize=false,
                        fsize=10,
                        qp1 = 0.05,
                        qp2 = 0.95,
                        labels= nothing,
                        doplot=true)
    hists, CI = makenuisancehists(optn, qp1, qp2, burninfrac = burninfrac,
                                     nbins = nbins, temperaturenum = temperaturenum,
                                 )
    if doplot
        fig,ax = subplots(length(hists), 1, figsize=figsize, squeeze=false)
        for (i, h) = enumerate(hists)
            bwidth, bx, denom = getbinsfromhist(h, pdfnormalize=pdfnormalize)
            ax[i].bar(bx, h.weights./denom, width=bwidth, edgecolor="black")
            !isnothing(labels) && ax[i].set_xlabel(labels[i])
            ax[i].set_ylabel("pdf")
        end
        nicenup(fig, fsize=fsize)
    end
    hists, CI
end

function getbinsfromhist(h, ;pdfnormalize=false)
    bwidth = diff(h.edges[1])
    bx = h.edges[1][1:end-1] + bwidth/2
    denom = pdfnormalize ? maximum(h.weights) : sum(h.weights)*bwidth
    bwidth, bx, denom
end        

function makenuisancehists(optn::OptionsNuisance, qp1, qp2; burninfrac = 0.5, nbins = 50,
    temperaturenum = -1)
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

function trapezoidrule(f,t)
    dt = diff(t)
    g  = 0.5*(f[1:end-1]+f[2:end])
    # integrated, tstart, tmid, tend
    # I recommend using tend
    cumsum(g.*dt), t[1:end-1], 0.5*(t[1:end-1]+t[2:end]), t[2:end]
end

function dodensityestimate(usekde::Bool, data, K::KDEtype, edges)
    if usekde
        kdefunc = kde_(K, data)
        kdefunc(0.5(edges[2:end]+edges[1:end-1]))
    else
        w = fit(Histogram, data, edges).weights
        w/(sum(w)*diff(edges)[1])
    end    
end

function stretchexists(F::Operator)
    in(:stretch, fieldnames(typeof(F))) && F.stretch
end

kde_(K::SJ, data) = kde_sj(;data)
kde_(K::LSCV, data) = kde_cv(;data)

abstract type KDEstimator end

# Sheather-Jones plugin
struct kde_sj <: KDEstimator
    data
    bw :: Real
end    

kde_sj(;data=zeros(0)) = kde_sj(data, bwsj(data))

function (foo::kde_sj)(xvals)
    density(foo.data, foo.bw, xvals)
end

function kde_sj(x; npoints=40)
    points = range(extrema(x)...,length=npoints)
    datapdf = kde_sj(data=x)
    datapdf(points), points
end

# LSCV KDE
struct kde_cv <: KDEstimator
    U
end    

kde_cv(;data=zeros(0)) = kde_cv(kde_lscv(data))

function (foo::kde_cv)(xvals)
    pdf(foo.U, xvals)
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
                usekde = false,
                kdetype = SJ(),
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
            himage[ilayer,:] = dodensityestimate(usekde, mthislayer, kdetype, edges)
            CI[ilayer,:] = [quantile(mthislayer,(qp1, 0.5, qp2))...]
            meanimage[ilayer] = mean(vec(mthislayer))
        else # there is an affine stretch with depth
            mthislayer = [m[ilayer] for m in M]
            expandedthislayer = T.(F.low[ilayer] .+ mthislayer*F.Δ[ilayer])
            himage[ilayer,:] = dodensityestimate(usekde, expandedthislayer, kdetype, edges)
            CI[ilayer,:] = [quantile(expandedthislayer,(qp1, 0.5, qp2))...]
            meanimage[ilayer] = mean(vec(expandedthislayer))
        end        
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


# some plotting codes for  AEM
function stepmodel(ax, iaxis, color, ρ , aem, model_lw, alpha)
    nfixed = aem.nfixed
    if isnothing(color) 
        ax[iaxis].step(ρ, aem.z[nfixed+1:end], linewidth=model_lw, alpha=alpha)
    else
        ax[iaxis].step(ρ, aem.z[nfixed+1:end]; linewidth=model_lw, alpha=alpha, color)
    end    
end   

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

function getsoundingsperline(soundings::Array{S, 1}) where S<:Sounding
    linestartidx = splitsoundingsbyline(soundings)                    
    nlines = length(linestartidx)                   
    s = map(1:nlines) do i
        a, b = linestartend(linestartidx, i, nlines, soundings)
        soundings[a:b]
    end    
end    

function linestartend(linestartidx, i, nlines, soundings)
    a = linestartidx[i]
    b = i != nlines ?  linestartidx[i+1]-1 : length(soundings)
    a, b
end    

function compatidxwarn(idx, lnames) 
    if !isempty(idx)
        if typeof(idx) == Array{Int64, 1} # if old format
            @warn "idx is same across ALL lines, not specific to line"
        else
            @assert !isnothing(lnames)
            @assert typeof(idx)==Vector{Vector{Int64}} "must be array of arrays per line"
        end    
    end
end

function docontinue(lnames, idx, soundings, a, b)
    continueflag = false
    idspec = []
    if !isempty(lnames) # only specific lines wanted, empty means all lines
        doesmatch = findfirst(lnames .== soundings[a].linenum) 
        if isnothing(doesmatch) 
            continueflag = true # do continue
        else
            @info lnames[doesmatch]
            if !isempty(idx)
                @assert length(lnames) == length(idx)
                @show idspec = idx[doesmatch]
                for id in idspec
                    @info "X, Y = $(soundings[a:b][id].X), $(soundings[a:b][id].Y)"
                end
            end    
        end    
    else # idx for the entire array of soundings
        idspec = idx
    end
    return continueflag, idspec 
end

function writevtkfromsounding(lineofsoundings::Array{S, 1}, zall) where S<:Sounding
    X, Y, Z = map(x->getfield.(lineofsoundings, x), (:X, :Y, :Z))
    lnum = lineofsoundings[1].linenum
    @info("opening summary: Line $(lnum)")
    rholow, rhomid, rhohigh = map(x->readdlm(x*"_line_$(lnum)_"*"summary.txt"), 
                                    ["rho_low", "rho_mid", "rho_hi"])
    writevtkfromxyzrho(rholow, rhomid, rhohigh, X, Y, Z, zall, lnum)
end

function getprobabilisticlinesfromdirectory(src_dir)
    fn = readdir(glob"rho_low*_summary_xyzrho.txt", src_dir)
    [parse(Int, split(split(f, "_summary")[1], "line_")[2]) for f in fn]
end

function writevtkfromxyzrho(rholow, rhomid, rhohigh, X, Y, Z, zall, lnum; dst_dir="")
    Ni, Nj = map(x->length(x), (X, Y))
    Nk = length(zall)
    x = [X[i] for i = 1:Ni, j = 1:1, k = 1:Nk]
    y = [Y[i] for i = 1:Ni, j = 1:1, k = 1:Nk]
    z = [Z[i] - zall[k] for i = 1:Ni, j = 1:1, k = 1:Nk]
    σlow, σmid, σhigh = map((rhohigh, rhomid, rholow)) do rho
        # switch from rho to sigma so low, hi interchanged 
        [-rho[k, i] for i = 1:Ni, j = 1:1, k = 1:Nk]
    end
    vtk_grid(joinpath(dst_dir, "Line_$(lnum)"), x, y, z) do vtk
        vtk["cond_low"]  = σlow
        vtk["cond_mid"]  = σmid
        vtk["cond_high"] = σhigh
    end
    nothing
end    

function writevtkphifromsummary(phid, sdev_phid, X, Y, Z, lnum; dst_dir="")
    Ni, Nj = map(x->length(x), (X, Y))
    x = [X[i] for i = 1:Ni, j = 1:1, k = 1:1]
    y = [Y[i] for i = 1:Ni, j = 1:1, k = 1:1]
    z = [Z[i] for i = 1:Ni, j = 1:1, k = 1:1]
    ϕmean, ϕsdev = map((phid, sdev_phid)) do p
        [p[i] for i = 1:Ni, j = 1:1, k = 1:1]
    end
    vtk_grid(joinpath(dst_dir, "phid_Line_$(lnum)"), x, y, z) do vtk
        vtk["phid_mean"]  = ϕmean
        vtk["phid_sdev"]  = ϕsdev
    end
    nothing
end    

function writevtkfromsounding(s::Vector{Array{S, 1}}, zall) where S<:Sounding
    pmap(s) do x
        writevtkfromsounding(x, zall)
    end    
end    

function readcols(cols::Vector, fname::String; decfactor=1, startfrom=1, dotill=nothing)
    d = readlargetextmatrix(fname, startfrom, decfactor, dotill)
    map(cols) do n
        # take signs of column numbers into account
        if (isa(n, Array))
            sign(n[1])*d[:,abs(n[1]):abs(n[2])]
        else
            if isa(n, Integer)
                sign(n)*d[:,abs(n)]
            elseif isa(n, Real)
                # pass through value
                fill(n, size(d,1))
            else
                @error "unknown entry type"  
            end        
        end    
    end    
end    

function getvtkwholeextent(fname)
    # assuming structured grid
    vtk = VTKFile(fname)
    dataset_element = ReadVTK.LightXML.root(vtk.xml_file)[vtk.file_type][1]
    whole_extent = parse.(Int, split(ReadVTK.LightXML.attribute(dataset_element, "WholeExtent", required = true), ' '))
    whole_extent
end    

function writevtkxmlforcurtain(vtkfname::String; src_epsg=0, dst_epsg=0, suffix="", vmin=0, vmax=0)
    whole_extent = getvtkwholeextent(vtkfname)
    nrows = whole_extent[2]
    ncols = whole_extent[6]
    writevtkxmlforcurtain(;vtkfname, nrows, ncols, src_epsg, dst_epsg, suffix, vmin, vmax)
end

function writevtkxmlforcurtain(;vtkfname="", 
    nrows=0, ncols=0, src_epsg=0, dst_epsg=0, suffix="", vmin=0, vmax=0)
    # nrows are usually ndepths
    # ncols are usually nsoundings
    # 0 based indexing so -1 than actual
    out =
    """
    <Layer version="1" layerType="VtkLayer">
        <DisplayName />
        <dimensions>0 $(nrows) 0 0 0 $(ncols)</dimensions>
        <URL>$(basename(vtkfname)*suffix)</URL>
        <sourceProjection>EPSG:$src_epsg</sourceProjection>
        <targetProjection>EPSG:$dst_epsg</targetProjection>
        <dataReader>vtkXmlReader</dataReader>
        <DefaultColourMap>Turbo</DefaultColourMap>
        <groupEvents>true</groupEvents>
        <minMaxOverride>true</minMaxOverride>
        <minOverride>$vmin</minOverride>
        <maxOverride>$vmax</maxOverride>
        <unitsOfMeasure>Log10 S/m</unitsOfMeasure>
    </Layer>
    """
    outname = split(vtkfname,".vts")[1]*".xml"
    f = open(outname, "w")
    write(f, out)
    close(f)
end

function colstovtk(cols::Vector, fname::String; decfactor=1, hasthick=true, islog10=false, prefix="")
    # for one file with one line, not used much I don't think 
    # TBD: delete function
    X, Y, Z, σ, thick = readcols(cols, fname; decfactor)
    thick = thick[1,:]
    zall = thicktodepth(thick; hasthick)
    fstring = basename(fname)
    dstring = dirname(fname)
    fpath = joinpath(dstring, "LEI_Line_"*prefix*fstring)
    colstovtk(X, Y, Z, σ, zall, fpath, islog10)
end

function colstovtk(X, Y, Z, σ, zall, fpath::String, islog10::Bool)
    # lowest level call
    Ni = length(X)
    Nk = length(zall)
    x = [X[i] for i = 1:Ni, j = 1:1, k = 1:Nk]
    y = [Y[i] for i = 1:Ni, j = 1:1, k = 1:Nk]
    z = [Z[i] - zall[k] for i = 1:Ni, j = 1:1, k = 1:Nk]
    σvtk = [σ[i, k] for i = 1:Ni, j = 1:1, k = 1:Nk]
    T = islog10 ?  x->x : x->log10(x) # if log 10 leave alone
    vtk_grid(fpath, x, y, z) do vtk
        vtk["cond_LEI"]  = T.(σvtk)
    end
    nothing
end

function colstovtk(cols::Dict, fname::String; decfactor=1, hasthick=true, islog10=false, prefix="")
    # dictionary call for a file with multiple lines
    Xc, Yc, Zc, σc, thickc, linec = map(x->get(cols, x, 0), ["X", "Y", "Z", "cond", "thick", "line"])
    X, Y, Z, σ, thick, lines = readcols([Xc, Yc, Zc, σc, thickc, linec], fname; decfactor)
    colstovtk(X, Y, Z, σ, thick, lines, fname; hasthick, islog10, prefix)
end    

function colstovtk(X, Y, Z, σ, thick_in, lines, fname::String; hasthick=true, islog10=false, prefix="")
    # no dictionary version next level call
    linenos  = unique(Int.(lines))
    if ndims(thick_in) == 2 
        thick = thick_in[1,:] # assumes all thicknesses are same
    else
        thick = thick_in
    end
    zall = thicktodepth(thick; hasthick)
    dstring = dirname(fname)
    for l in linenos
        idx = lines .== l
        fstring = "$l"
        @info "doing Line "*fstring
        fpath = joinpath(dstring, "LEI_Line_"*prefix*fstring)
        colstovtk(X[idx], Y[idx], Z[idx], σ[idx,:], zall, fpath, islog10)
    end
end

function thicktodepth(thick; hasthick=true)
    if hasthick
        zb = [0; cumsum(thick)]
        zall = 0.5(zb[1:end-1]+zb[2:end])
    else
        zall = thick    
    end
    zall    
end

function zcentertoboundary(zall)
    zb = zeros(length(zall))
    zb[1] = 0.
    zb[2] = 2*zall[1]
    for i in 3:length(zb)
        delz = 2(zall[i-1] - zb[i-1])
        zb[i] = zb[i-1] + delz
    end    
    zb
end

function zboundarytocenter(zb)
    # only for exact depth of last layer halfspace for geometric thickness increase
    # default used in transD_GP. If you have other generic thickness
    # use thicktodepth()
    
    # first get extendfrac r
    numerator = diff(zb[2:end])
    denom = diff(zb[1:end-1])
    r = (denom'*denom)\(denom'*numerator) # overkill but I love least squares
    # now for dz
    numerator = diff(zb)
    denom = map(2:length(zb)) do n
        r^(n-2)
    end
    dz = (denom'*denom)\(denom'*numerator) # also overkill
    zall,  = setupz(0.0, r; dz, n=length(zb))
    zall
end  

function writeijkfromsounding(s::Vector{Array{S, 1}}, zall) where S<:Sounding
    pmap(s) do x
        writeijkfromsounding(x, zall)
    end    
end 

function writeijkfromsounding(lineofsoundings::Array{S, 1}, zall) where S<:Sounding
    X, Y, Z = map(x->getfield.(lineofsoundings, x), (:X, :Y, :Z))
    lnum = lineofsoundings[1].linenum
    @info("opening summary: Line $(lnum)")
    rholow, rhomid, rhohigh = map(x->readdlm(x*"_line_$(lnum)_"*"summary.txt"), 
                                    ["rho_low", "rho_mid", "rho_hi"])
    Ni, Nj = map(x->length(x), (X, Y))
    Nk = length(zall)
    x = [X[i] for i = 1:Ni, j = 1:1, k = 1:Nk]
    y = [Y[i] for i = 1:Ni, j = 1:1, k = 1:Nk]
    z = [Z[i] - zall[k] for i = 1:Ni, j = 1:1, k = 1:Nk]
    i = [i-1 for i = 1:Ni, j = 1:1, k = 1:Nk]
    j = [k-1 for i = 1:Ni, j = 1:1, k = 1:Nk]
    k = zeros(size(j))
    σlow, σmid, σhigh = map((rhohigh, rhomid, rholow)) do rho
        # switch from rho to sigma so low, hi interchanged 
        [-rho[k, i] for i = 1:Ni, j = 1:1, k = 1:Nk]
    end
    str = ["low", "mid", "high"]
    map(zip(str, [σlow, σmid, σhigh ])) do (f, σ)
        writeijkfromgrid(f, lnum, σ, x, y, z, i, j, k, Ni, Nk)
    end    
    nothing
end

function writeijkfromgrid(str, lnum, σ, x, y, z, i, j, k, Ni, Nk)
    fname_data = "$(lnum)_"*str*".sg.data"
    fname_header = fname_data[1:end-5]
    # write_data
    open(fname_data, "w") do f
        write(f, "*\n")
        write(f, "*   X   Y   Z  log10Spm  I   J   K\n")
        write(f, "*\n")
    end
    io = open(fname_data, "a")
    for c in 1:length(σ)
        msg = @sprintf("%.2f %.2f %.2f %.4f %10i %10i %10i\n", x[c], y[c], z[c], σ[c], i[c], j[c], k[c])
        write(io, msg)
    end
    close(io)
    # write header
    io = open(fname_header, "w")
    headerlines = """
    GOCAD SGrid 1
    HEADER {
    name:$(lnum)
    painted:true
    *painted*variable:log10Cond
    ascii:on
    double_precision_binary:off
    cage:false
    volume:true
    *volume*grid:false
    *volume*transparency_allowed:false
    *volume*points:false
    shaded_painted:false
    precise_painted:true
    *psections*grid:false
    *psections*solid:true
    dead_cells_faces:false
    *metadata*Line: $lnum
    *metadata*Type: $str
    }

    AXIS_N $(Ni) $(Nk) 1
    PROP_ALIGNMENT POINTS
    ASCII_DATA_FILE $(fname_data)

    PROPERTY 1 log10Cond
    PROP_UNIT 1 log10Spm
    PROP_NO_DATA_VALUE 1 -999

    END
    """
    write(io, headerlines)
    close(io)
    nothing    
end

function findclosestidxincolfile(Xwanted, Ywanted, cols::Vector, fname::String; decfactor=1, hasthick=true)
    X, Y, σ, thick = readcols(cols, fname; decfactor)
    thick = thick[1,:]
    zall = thicktodepth(thick; hasthick)
    zb = zcentertoboundary(zall)
    XY = [X';Y']
    idx = getclosestidx(XY, Xwanted, Ywanted)
    σ[idx,:], zb
end

function makegrid(vals::AbstractArray, soundings::Array{S, 1}; donn=false,
    dr=10, zall=[NaN], dz=-1) where S<:Sounding
    X = [s.X for s in soundings]
    Y = [s.Y for s in soundings]
    topo = [s.Z for s in soundings]
    makegrid(vals, X, Y, topo; donn, dr, zall, dz)
end

function makegrid(vals::AbstractArray, X, Y, topo; donn=false,
    dr=10, zall=[NaN], dz=-1)
    @assert all(.!isnan.(zall)) 
    @assert dz>0
    R, gridr = getRandgridr(X, Y, dr)
    # height = topo in mAHD - depth # mAHD
    minahd = minimum(topo) - maximum(zall) # should really be zboundary[end] but are both in halfspace
    maxahd = maximum(topo) # top is the one to get right
    gridz = range(maxahd, minahd, step=-dz) 
    gridr_, gridz_ = map((gridr, gridz)) do a 
        0.5(a[1:end-1]+a[2:end])
    end
    if donn
        rr, zz = [r for z in zall, r in R], [z for z in zall, r in R]
        zz = topo' .- zz # mAHD
        kdtree = KDTree([rr[:]'; zz[:]'])
        rr, zz = [r for z in gridz_, r in gridr_], [z for z in gridz_, r in gridr_]
        idxs, = nn(kdtree, [rr[:]'; zz[:]'])
        img = zeros(size(rr))
        for i = 1:length(img)
            img[i] = vals[idxs[i]]
        end
        topofine = gridpoints(R, gridr_, topo)
        topoout = gridpoints(R, gridr, topo)
    else
        nodes = ([z for z in zall], [r for r in R])
        itp = extrapolate(interpolate(nodes, vals, Gridded(Linear())), Line()) 
        topofine = (interpolate((R,), topo, Gridded(Linear())))(gridr_)
        topoout = (interpolate((R,), topo, Gridded(Linear())))(gridr)
        img = [itp(topofine[iy] - x,y) for x in gridz_, (iy, y) in enumerate(gridr_)]
        zz = [z for z in gridz_, r in gridr_] 
    end    
    img[zz .>topofine'] .= NaN
    img[zz .< topofine' .- maximum(zall)] .= NaN # should be zboundary[end] but both in halfspace
    img, gridr, gridz, topoout, R
end

function getRandgridr(X, Y, dr; rangelenth=0)
    R = cumulativelinedist(X, Y)
    gridr = rangelenth == 0 ? range(R[1], R[end], step=dr) : range(R[1], R[end], length=rangelenth)
    R, gridr
end

function getnextxyinr(XY1, XY2, Δr)
    Δx, Δy = XY2 - XY1
    r̂  = normalize([Δx, Δy])
    @assert sqrt(Δx^2 + Δy^2) >= Δr "Δr throws it beyond next point"
    XY1 + Δr*r̂
end

function getXYlast(X, Y, Δr; rangelenth=0)
    R, gridredges = getRandgridr(X, Y, Δr; rangelenth)
    enddist = R[end] - gridredges[end]
    Xend, Yend = X[end], Y[end]
    if enddist > 0 
        remdist = gridredges[end]-R[end-1]
        Xend, Yend = getnextxyinr([X[end-1], Y[end-1]], [X[end], Y[end]], remdist)
    end    
    R[end] = gridredges[end]
    Xexact = vcat(X[1:end-1], Xend)
    Yexact = vcat(Y[1:end-1], Yend)
    Xexact, Yexact, R, gridredges
end

function getallxyinr(Xin, Yin, Δr; rangelenth=0)
    # the above function returns middle of points in gridr
    X, Y, R, gridredges = getXYlast(Xin, Yin, Δr; rangelenth)
    # these are cell centres
    gridr = 0.5(gridredges[1:end-1]+gridredges[2:end])
    xyfine = zeros(2, length(gridr))
    c = 1
    xyfine[:,1] = getnextxyinr([X[1], Y[1]], [X[c+1], Y[c+1]], Δr/2)
    for i = 1:length(gridr)-1
        if gridr[i+1]<R[c+1]
            xyfine[:,i+1] = getnextxyinr(xyfine[:,i], [X[c+1], Y[c+1]], Δr)
        else
            Δrsmall = gridr[i+1]-R[c+1]
            xyfine[:,i+1] = getnextxyinr([X[c+1], Y[c+1]], [X[c+2], Y[c+2]], Δrsmall)
            c += 1
        end
    end
    xyfine, gridr
end
    
function cumulativelinedist(X,Y)
    dx = diff(X)
    dy = diff(Y)
    R = [0.; cumsum(sqrt.(dx.^2 + dy.^2))]
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
    ΔX = X[end] - X[1]
    ΔY = Y[end] - Y[1]
    NSline, EWline = false, false
    if abs(ΔX/ΔY) < 0.05
        NSline = true
    elseif abs(ΔY/ΔX) < 0.05
        EWline = true
    end    
    
    Eislast, Nislast = true, true
    ΔX<0 && (Eislast = false)
    ΔY<0 && (Nislast = false)

    Eislast, Nislast, EWline, NSline
end

function makesummarygrid(soundings, pl, pm, ph, ρmean, zall, dz; dr=10)
    # first flip ρ to σ and so pl and ph are interchanged
    phgrid, gridx, gridz, topofine, R = makegrid(-pl, soundings, zall=zall, dz=dz, dr=dr)
    plgrid,                           = makegrid(-ph, soundings, zall=zall, dz=dz, dr=dr)
    pmgrid,                           = makegrid(-pm, soundings, zall=zall, dz=dz, dr=dr)
    σmeangrid,                        = makegrid(-ρmean, soundings, zall=zall, dz=dz, dr=dr)
    Z = [s.Z for s in soundings]
    phgrid, plgrid, pmgrid, σmeangrid, gridx, gridz, topofine, R, Z
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

function plotNEWSlabels(soundings, gridx, gridz, axarray, 
        x0=nothing, y0=nothing, xend=nothing, yend=nothing;
        preferEright=false, preferNright=false, fontsize=12)
    Eislast, Nislast, EWline, NSline = whichislast(soundings)
    beginpos, endpos = "", ""
    if !any(isnothing.([x0,y0,xend,yend]))
        beginpos = @sprintf(" %.1f", x0/1000)*@sprintf(" %.1f", y0/1000)
        endpos = @sprintf(" %.1f", xend/1000)*@sprintf(" %.1f", yend/1000)
    end
    for s in axarray
        minylim = minimum(s.get_ylim())
        inverted = false
        if (preferNright && !Nislast) || (preferEright && !Eislast) 
            inverted = true
        end 
        
        if Eislast & Nislast # old WN
            dir1, dir2 = "SW", "NE"
            if EWline
                dir1, dir2 = "W", "E"
            elseif NSline
                dir1, dir2 = "S", "N"
            end    
        elseif !Eislast & Nislast # old EN
            dir1, dir2 = "SE", "NW"
            if EWline
                dir1, dir2 = "E", "W"
            elseif NSline
                dir1, dir2 = "S", "N"
            end 
        elseif Eislast & !Nislast # old WS
            dir1, dir2 = "NW", "SE"
            if EWline
                dir1, dir2 = "W", "E"
            elseif NSline
                dir1, dir2 = "N", "S"
            end 
        else # old ES
            dir1, dir2 = "NE", "SW"
            if EWline
                dir1, dir2 = "E", "W"
            elseif NSline
                dir1, dir2 = "N", "S"
            end 
        end
        ha = inverted ? "right" : "left"
        s.text(gridx[1], minylim, dir1*beginpos, backgroundcolor=s.get_facecolor(), ha = ha; fontsize)
        ha = inverted ? "left" : "right"
        s.text(gridx[end], minylim, dir2*endpos, backgroundcolor=s.get_facecolor(), ha = ha; fontsize)      
    end
    flipline(preferNright, preferEright, Nislast, Eislast, axarray[end])
end

function flipline(preferNright, preferEright, Nislast, Eislast, ax)
    @assert !(preferEright & preferNright) #both true means don't know what you're doing
    (preferNright && !Nislast) && ax.invert_xaxis()
    (preferEright && !Eislast) && ax.invert_xaxis()
end

function plotprofile(ax, idxs, Z, R)
    for idx in idxs
        ax.plot(R[idx]*[1,1], [ax.get_ylim()[1], Z[idx]], "-w")
        ax.plot(R[idx]*[1,1], [ax.get_ylim()[1], Z[idx]], "--k")
    end
end

function getRsplits(R, Rmax)
    q, rmn = divrem(maximum(R), Rmax)
    thereisrmn = !iszero(rmn)
    idx = thereisrmn ? zeros(Int, Int(q+1)) : zeros(Int, Int(q)) 
    i = 1
    for (ir, r) in enumerate(R)
        if r >= Rmax*i
            idx[i] = ir
            i += 1
        end
    end
    idx, thereisrmn        
end    

function getinterpsplits(idx_split, nimages, gridx)
    i_idx = 1:nimages
    ab = zeros(Int, nimages,2)
    for i in i_idx
        a = i == firstindex(i_idx) ? 1 : idx_split[i-1]
        b = i != lastindex(i_idx)  ? idx_split[i] : lastindex(gridx)
        b = b-1 # as we are now providing 1 less than the number of edges since 26c19a0d2e4a
        ab[i,1] = a
        ab[i,2] = b
    end
    if ab[end,1] == ab[end,2]
        ab = ab[1:end-1,:]
        nimages = nimages-1
    end
    ab, nimages
end

function plotsummarygrids1(soundings, meangrid, phgrid, plgrid, pmgrid, gridx, gridz, topofine, R, Z, χ²mean, χ²sd, lname; qp1=0.05, qp2=0.95,
                        figsize=(10,10), fontsize=12, cmap="turbo", vmin=-2, vmax=0.5, Rmax=nothing,
                        topowidth=2, idx=nothing, omitconvergence=false, useML=false, preferEright=false, preferNright=false,
                        saveplot=false, yl=nothing, dpi=300, showplot=true, showmean=false, logscale=true)
    if isnothing(Rmax)
        Rmax = maximum(gridx)
    end
    while true
        idx_split, thereisrmn = getRsplits(gridx, Rmax)
        nimages = length(idx_split)
        dr = gridx[2] - gridx[1]
        dr >= Rmax && error("dr:$dr >= Rmax:$Rmax, decrease dr")
        if iszero(idx_split[1]) && thereisrmn # Rmax is larger than section
            nx = length(range(gridx[1], Rmax, step=dr))
        elseif !iszero(idx_split[1]) && !thereisrmn # Rmax is exactly the section length
            nx = length(gridx)
        elseif iszero(idx_split[2]) # Rmax is smaller than the section
            nx = idx_split[1]
        else    
            nx = idx_split[2]-idx_split[1] # There are many Rmax length splits 
        end    
        ab, nimages = getinterpsplits(idx_split, nimages, gridx)
        i_idx = 1:nimages
        for i in i_idx
            a, b = ab[i,1], ab[i,2]
            a_uninterp = i == firstindex(i_idx) ? 1 : findlast(R.<=gridx[a])
            b_uninterp = i != lastindex(i_idx)  ? findlast(R.<=gridx[b]) : lastindex(soundings)
            if thereisrmn && i == lastindex(i_idx)
                xrangelast = range(gridx[a], step=dr, length=nx)
            else 
                xrangelast = nothing
            end
            f, s, icol = setupconductivityplot(gridx[a:b], omitconvergence, showmean, R[a_uninterp:b_uninterp], 
                figsize, fontsize, lname, χ²mean[a_uninterp:b_uninterp], χ²sd[a_uninterp:b_uninterp], useML, i, nimages, logscale)
            
            summaryconductivity(s, icol, f, soundings[a_uninterp:b_uninterp], 
                meangrid[:,a:b], phgrid[:,a:b], plgrid[:,a:b], pmgrid[:,a:b], 
                gridx[a:b], gridz, topofine[a:b], R[a_uninterp:b_uninterp], Z[a_uninterp:b_uninterp], ; qp1, qp2, fontsize, 
                cmap, vmin, vmax, topowidth, idx, omitconvergence, preferEright, preferNright, yl, showmean, xrangelast)
            postfix = nimages == 1 ? "_whole" : "_split_$(i)_of_$(nimages)"
            saveplot && savefig(lname*postfix*".png", dpi=dpi)
            showplot || close(f)
        end
        nimages == 1 && break # don't need another pass if it was whole line to start with
        Rmax = maximum(gridx) # this will ensure nimages=1 on next pass
    end
end

function setupconductivityplot(gridx, omitconvergence, showmean, R, figsize, fontsize, lname, χ²mean, χ²sd, useML, iimage, nimages, logscale)
    dr = diff(gridx)[1]
    nrows = omitconvergence ? 5 : 6
    height_ratios = omitconvergence ? [1, 1, 1, 1, 0.1] : [0.4, 1, 1, 1, 1, 0.1]
    if !showmean 
        height_ratios = [height_ratios[1:2]..., height_ratios[4:end]...]
        nrows-=1
    end    
    f, s = plt.subplots(nrows, 1, gridspec_kw=Dict("height_ratios" => height_ratios),
                        figsize=figsize)
    f.suptitle(lname*" Δx=$dr m, Fids: $(length(R)), $(iimage) of $(nimages)", fontsize=fontsize)
    icol = 1
    if !omitconvergence
        s[icol].plot(R, χ²mean)
        linecolor = [s[icol].get_facecolor()...]
        linecolor[1:3] = abs.(linecolor[1:3] .- 1)
        s[icol].plot(R, ones(length(R)), "--", color=linecolor)
        s[icol].fill_between(R, vec(χ²mean-χ²sd), vec(χ²mean+χ²sd), alpha=0.5)
        s[icol].set_ylabel(L"ϕ_d")
        titlestring = useML ? "Max likelihood variance adjustment" : "Data misfit"
        s[icol].set_title(titlestring)
        logscale && s[icol].set_yscale("log")
        icol += 1
    end
    f, s, icol
end    

function summaryconductivity(s, icol, f, soundings, meangrid, phgrid, plgrid, pmgrid, gridx, gridz, topofine, R, Z, ; qp1=0.05, qp2=0.95,
    fontsize=12, cmap="turbo", vmin=-2, vmax=0.5, topowidth=2, idx=nothing, omitconvergence=false, preferEright=false, preferNright=false,
    yl=nothing, showmean=false, xrangelast=nothing)
    icolstart = icol
    s[icol].imshow(plgrid, cmap=cmap, aspect="auto", vmax=vmax, vmin = vmin,
                extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
    s[icol].plot(gridx, topofine, linewidth=topowidth, "-k")
    idx == nothing || plotprofile(s[icol], idx, Z, R)
    s[icol].set_title("Percentile $(round(Int, 100*qp1)) conductivity")
    s[icol].set_ylabel("Height m")
    omitconvergence || s[icol].sharex(s[icol-1])
    icol += 1
    s[icol].imshow(pmgrid, cmap=cmap, aspect="auto", vmax=vmax, vmin = vmin,
                extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
    s[icol].plot(gridx, topofine, linewidth=topowidth, "-k")
    s[icol].set_title("Percentile 50 conductivity")
    idx == nothing || plotprofile(s[icol], idx, Z, R)
    s[icol].set_ylabel("Height m")
    s[icol].sharex(s[icol-1])
    s[icol].sharey(s[icol-1])
    icol += 1
    if showmean
        s[icol].imshow(meangrid, cmap=cmap, aspect="auto", vmax=vmax, vmin = vmin,
                    extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
        s[icol].plot(gridx, topofine, linewidth=topowidth, "-k")
        s[icol].set_title("Mean conductivity")
        s[icol].set_ylabel("Height m")
        idx == nothing || plotprofile(s[icol], idx, Z, R)
        s[icol].sharex(s[icol-1])
        s[icol].sharey(s[icol-1])
        icol +=1
    end    
    imlast = s[icol].imshow(phgrid, cmap=cmap, aspect="auto", vmax=vmax, vmin = vmin,
                extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
    s[icol].plot(gridx, topofine, linewidth=topowidth, "-k")
    s[icol].set_xlabel("Line distance m", labelpad=0)
    s[icol].set_title("Percentile $(round(Int, 100*qp2)) conductivity")
    s[icol].set_ylabel("Height m")
    idx == nothing || plotprofile(s[icol], idx, Z, R)
    s[icol].sharex(s[icol-1])
    s[icol].sharey(s[icol-1])
    s[icol].set_xlim(extrema(gridx))
    # map(x->x.set_xticklabels([]), s[1:end-2])
    map(x->x.tick_params(labelbottom=false), s[1:end-2])
    # map(x->x.grid(), s[1:end-1])
    isa(yl, Nothing) || s[end-1].set_ylim(yl...)
    x0, y0 = soundings[1].X, soundings[1].Y
    xend, yend = soundings[end].X, soundings[end].Y
    if !isnothing(xrangelast)
        s[icol].set_xlim(extrema(xrangelast))
    end    
    plotNEWSlabels(soundings, gridx, gridz, s[icolstart:end-1], x0, y0, xend, yend; preferEright, preferNright, fontsize)
    cb = f.colorbar(imlast, cax=s[end], orientation="horizontal")
    cb.set_label("Log₁₀ S/m", labelpad=0)
    nicenup(f, fsize=fontsize, h_pad=0)
    label = f._suptitle.get_text()
    VE = round(Int, getVE(s[end-1]))
    f.suptitle(label*", VE=$(VE)X", fontsize=0.8*fontsize)
end

function getVE(ax)
    figW, figH = ax.get_figure().get_size_inches()
    # Axis size on figure
    _, _, w, h = ax.get_position().bounds
    # Ratio of display units
    disp_ratio = (figH * h) / (figW * w)
    # Ratio of data units
    # Negative over negative because of the order of subtraction
    data_ratio = abs.(diff([ax.get_ylim()...])[1] / diff([ax.get_xlim()...])[1])
    return disp_ratio / data_ratio
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
            Mblock[i,j] = sum(idxdepth) == 0 ? NaN : eval(cond)(m[idxdepth])
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

function getclosestidx(Xwell, Ywell, soundings::Vector{S}) where S<: Sounding
    XY = hcat([[s.X, s.Y] for s in soundings]...)
    getclosestidx(XY, Xwell, Ywell)
end    

function getclosestidx(XY, Xwanted, Ywanted; showinfo=true)
    idx, dist = getclosestidxanddist(XY, Xwanted, Ywanted)
    showinfo && @info "distance is $dist"
    idx
end    

function getclosestidxanddist(XY, Xwanted, Ywanted)
    tree = KDTree(XY)
    idx, dist = nn(tree, [Xwanted;Ywanted])
end

function calchighalt(fname_dat, windownum_start, windownum_end; 
                     times=nothing, figsize=(5,5), fontsize=10, time_units="", amp_units="",
                     suffix="", agreednoise=nothing, agreedbias=nothing)
    # reading and calculating high_alt
    if !isnothing(agreednoise)
        @assert (ndims(agreednoise) == 2) "we need envelope and times"
        @assert (ndims(agreedbias)  == 2) "we need envelope and times"
    end
    nwindows = windownum_end - windownum_start + 1
    if !isnothing(times)
        @assert length(times) = nwindows
        xlab = "time "*time_units
        dologx = true
    else
        times = 1:nwindows
        xlab = "windows"
        dologx = false
    end
    halt_data = readlargetextmatrix(fname_dat)[:,windownum_start:windownum_end]
    bias = vec(abs.(infnanmean(halt_data, 1)))
    sdev = vec(infnanstd(halt_data, 1))
    f, ax = plt.subplots(1, 1; figsize)
    ax.semilogy(times, bias, label="bias calculated")
    ax.semilogy(times, sdev, label="sdev calculated")
    if !isnothing(agreednoise)
        ax.semilogy(times, agreedbias, label="bias agreed")
        ax.semilogy(times, agreedsdev, label="sdev agreed")
    end
    ax.set_title(suffix*" noise")
    ax.set_xlabel(xlab)
    ax.set_ylabel(amp_units)
    dologx && ax.set_xcale("log")
    nicenup(f, fsize=fontsize)
    !isempty(suffix) && (suffix = "_"*suffix)
    savefig(fname_dat[1:end-4]*suffix*".png")
    writedlm(fname_dat[1:end-4]*suffix*".txt", sdev)
end

#function to read the *dfn file and extract the column number and column names as a *.txt file 
function dfn2hdr(dfnfile::String; writecorrectedhdr=false)
    dfn = readlines(dfnfile)
    dfn = dfn[(.!contains.(dfn, "RT=PROJ") .& .!contains.(dfn, "RT=TRNS")
              .& .!contains.(dfn, "RT=COMM"))]

    rectype_rgx = r"RT=([^;]*);"
    # name_rgx = r"NAME=" # we don't need the name field
    name_rgx = deepcopy(rectype_rgx)
    # inc_regex = r":([0-9]+)f"i #this will only work with floating point fields
    inc_regex = r":\s*([0-9]+)(f|e)"i # this works with spaces before colon and with exponential fields

    cumulative_columns = 0  #this will set up a cumulative variable 
    
    # the name of HDR file associated with DFN
    fname_hdr = getgdfprefix(dfnfile)*(writecorrectedhdr ? "_corrected.hdr" : ".hdr") 
    io = open(fname_hdr,"w")    
    data_first = false
    for row = dfn
        m = match(rectype_rgx, row)
        isnothing(m) && continue
        data_first |= (m[1] == "DATA")

        inc_match = match(inc_regex, row)  #matching the keyword with the string 
        if isnothing(inc_match) # single column field
            inc = 1
            idx2 = split(row, name_rgx)
            idx2_r = first.(split.(idx2[2], ":"))
            # @info idx2_r
        else # multiple column field
            inc = parse(Int64, inc_match.captures[1]) #the outcome is int64
            idx2 = split(row, rectype_rgx)
            idx2_r = first.(split.(idx2[2], ":"))
        end

        firstcol = cumulative_columns + 1 + data_first
        if writecorrectedhdr 
            firstcol -= 1
        end    
        lastcol = firstcol + inc - 1
        cumulative_columns += inc
        
        if occursin("END DEFN", idx2_r)  #type is string 
            idx2_r = first.(split.(idx2_r,";"))
            break # proper termination
        end

        #start writing to a file here
        if (firstcol != lastcol)
            writedlm(io,[firstcol lastcol idx2_r]) 
        else
            writedlm(io, [firstcol idx2_r])
        end
    end
    close(io)
end

function correctOMaseggdf(fname::String; writecorrectedhdr=false)
    # use this, to read ASEG-GDF generated by Oasis Montaj (upto at least 2023.2),
    # which has an incorrect first column of text. This function writes out the correct
    # .DAT file omitting the column and the corrected HDR file withe the correct headers 
    # and corresponding column numbers

    f = open(fname)
    fcorr = open(fname[1:end-4]*"_corrected.dat", "w")
    for (i, str) in enumerate(eachline(f))
        if startswith(str, "DATA")
            newline = replace(str, "DATA" => "", count=1)*"\n"
            write(fcorr, newline)
        end    
    end
    close.([f, fcorr])
    if writecorrectedhdr
        smallfile, capsfile = fname[1:end-3].*["DFN", "dfn"]
        if isfile(smallfile)
            file = smallfile
        elseif isfile(capsfile)
            file = capsfile
        else
            @error "file not found"
        end
        dfn2hdr(file; writecorrectedhdr)
    end    
    nothing
end    

function getgdfprefix(dfnfile::String)
    location = findfirst(".dfn", lowercase(dfnfile))[1] # get file prefix
    dfnfile[1:location-1] # the prefix
end    

function readlargetextmatrix(fname::String)
    a = read(fname) # as UInt8
    map!(c -> c == UInt8('\t') ? UInt8(' ') : c, a, a) # replace tab with space and use memory for a
    # convert DataFrame to Float64 matrix
    CSV.File(IOBuffer(a); ignorerepeated=true, types=Float64, header=false, delim=' ')|>CSV.Tables.matrix
end

function readlargetextmatrix(fname::String, startfrom, skipevery, dotillsounding::Union{Int, Nothing})
    soundings = readlargetextmatrix(fname)
    if !isnothing(dotillsounding)
        soundings = soundings[startfrom:skipevery:dotillsounding,:]
    else
        soundings = soundings[startfrom:skipevery:end,:]
    end
    soundings
end

function pairinteractionplot(d; varnames=nothing, figsize=(8.5,6), nbins=25, fontsize=8, fbounds=nothing,
        cmap="bone_r", islogpdf=false, showpdf=false, vecofpoints=nothing, vecofpointscolor=nothing,
        scattersize=1, scattercolor="yellowgreen", scatteralpha=1)
    # plot pairs of scatter, d assumed to have realisations in rows
    # vecofpoints should be a vector of vectors if nothing
    @assert !isnothing(varnames)
    nvars = size(d,2)
    f = figure(figsize=figsize)
    T = islogpdf ? x->log(x) : x->x
    for i=1:size(d,2)
        for j=1:i
            c = getrowwise(i,j,nvars)
            ax = subplot(nvars, nvars, c)
            if i==j 
                h = fit(Histogram, d[:,i], nbins=nbins)
                bwidth, bx, denom = getbinsfromhist(h, pdfnormalize=true)
                ax.bar(bx, h.weights./denom, width=bwidth, edgecolor="none", color="yellowgreen")
                ax.set_ylabel("probability")
                ax.yaxis.set_label_position("right")
                ax.tick_params(axis="y", labelright=true, labelleft=false, right=true, left=false)
                if !isnothing(vecofpoints)    
                    for (iv,v) in enumerate(vecofpoints)
                        if isnothing(vecofpointscolor) 
                            ax.plot(v[i]*ones(2),[0, 1], markeredgewidth=3)
                        else
                            ax.plot(v[i]*ones(2),[0, 1], color=vecofpointscolor[iv], linewidth=3)    
                        end
                    end
                end
                if !isnothing(fbounds)    
                    ax.plot(fbounds[i, 1]*ones(2),[0, 1], "--k", linewidth=1)
                    ax.plot(fbounds[i, 2]*ones(2),[0, 1], "--k", linewidth=1)   
                end            
            else
                h = fit(Histogram, (d[:,i], d[:,j]), nbins=nbins)
                if showpdf
                    ax.pcolormesh(h.edges[2], h.edges[1], T.(h.weights/maximum(h.weights)), cmap=cmap)
                else    
                    ax.scatter(d[:,j], d[:,i], s=scattersize, alpha=scatteralpha, c=scattercolor)
                end
                if !isnothing(vecofpoints)
                    for (iv,v) in enumerate(vecofpoints)
                        if isnothing(vecofpointscolor) 
                            ax.plot(v[j], v[i], "+", markersize=10*scattersize, markeredgewidth=3)
                        else
                            ax.plot(v[j], v[i], "+", color=vecofpointscolor[iv], 
                            markersize=10*scattersize, markeredgewidth=3)
                        end        
                    end
                end
                if !isnothing(fbounds)    
                    ax.plot(fbounds[j,1]*ones(2),fbounds[i,:], "--k", linewidth=1)
                    ax.plot(fbounds[j,2]*ones(2),fbounds[i,:], "--k", linewidth=1)
                    ax.plot(fbounds[j,:],fbounds[i,1]*ones(2), "--k", linewidth=1)
                    ax.plot(fbounds[j,:],fbounds[i,2]*ones(2), "--k", linewidth=1)
                end         
                j == 1 && ax.set_ylabel(varnames[i])
                j == 1 || ax.tick_params(labelleft=false)
            end
            (i == nvars) && ax.set_xlabel(varnames[j])
            ax.ticklabel_format(useOffset=false) 
            (i == nvars) ||  ax.tick_params(labelbottom=false)
        end
    end
    ax = f.axes
    for i = size(d,2)-1:-1:1
    # align x axes    
        for j = i:-1:1
            this = firstval(i) + j-1
            next = firstval(i+1) + j-1
            ax[this].sharex(ax[next])
        end
    end
    for i = size(d,2):-1:1
        # align y axes, abundance of caution for
        # weird histograms probably won't need this    
        for j = i-2:-1:1
            this = firstval(i) + j-1
            next = this + 1
            ax[this].sharey(ax[next])
        end
    end    
    # very last histogram shares x 
    # with y axis of previous
    ax[end].set_xlim(ax[end-1].get_ylim())
    nicenup(f, fsize=fontsize)
end

function plotgausshist(x;nbins=20, figsize=(4,4), title="")
    f, ax = plt.subplots(figsize=figsize)
    plotgausshist(x, ax;nbins)
    ax.set_title(title)
    nicenup(f)
end

function plotgausshist(x, ax;nbins=20, color=nothing)
    h = normalize(fit(Histogram, x; nbins), mode=:pdf)
    edges = h.edges[1]
    width=diff(edges)
    ax.bar(edges[1:end-1], h.weights, align="edge", width=width)
    if isnothing(color)
        ax.plot(edges, 1/sqrt(2pi)*exp.(-edges.^2), "--")
    else
        ax.plot(edges, 1/sqrt(2pi)*exp.(-edges.^2), "--", color=color)
    end
end

nanmean(x) = mean(filter(!isnan,x))
nanmean(x, dims) = mapslices(nanmean,x,dims=dims)
nanstd(x) = std(filter(!isnan,x))
nanstd(x, dims) = mapslices(nanstd, x, dims=dims)

infmean(x) = mean(filter(!isinf,x))
infmean(x, dims) = mapslices(infmean,x,dims=dims)
infstd(x) = std(filter(!isinf,x))
infstd(x, dims) = mapslices(infstd, x, dims=dims)

function infnanmean(x)
    idx = isnan.(x) .| isinf.(x)
    mean(x[.!idx])
end

infnanmean(x, dims) = mapslices(infnanmean,x,dims=dims)

function infnanstd(x)
    idx = .!isnan.(x) .& .!isinf.(x)
    std(x[idx])
end

infnanstd(x, dims) = mapslices(infnanstd,x,dims=dims)

getrowwise(i,j,nvars) = (i-1)*nvars+j
getcolwise(i,j,nvars) = (j-1)*nvars+i
firstval(n) = n == 1 ? 1 : firstval(n-1) + n-1 # index number when filling rowwise upto diagonals

# these are a little hacky for reading from gradient inversion and probabilistic files with minimal info
function readfzipped(fzipped::String, nlayers::Int; nnu=0)
    A = readdlm(fzipped)
    X, Y, Z, fid, line = map(i->A[:,i],(1:5))
    nu = A[:,end-nnu:end-1]
    zall = A[:,end-2nlayers-nnu:end-nlayers-nnu-1][1,:]
    σ = A[:,end-nlayers-nnu:end-nnu-1] # so we can plot TEMPEST and SPECTREM similar to heli
    ϕd = A[:,end]
    X, Y, Z, fid, line, zall, σ, ϕd, nu
end  

function readfzipped(fzipped::String, line::Int, nlayers::Int; nnu=0)
    X, Y, Z, fid, linesall, zall, σ, ϕd, nu = readfzipped(fzipped, nlayers; nnu)
    idx = linesall .== line
    @assert !isempty(idx)
    X, Y, Z, fid, linesall, σ, ϕd, nu = map(x->x[idx,:],(X, Y, Z, fid, linesall, σ, ϕd, nu))
    X, Y, Z, fid, linesall, zall, σ', ϕd, nu
end

function getdeterministicoutputs(outputs::AbstractArray) 
    # useful for plotting an aray of outputs read programmatically from readfzipped
    # same as 
    # X, Y, Z, fid, line, zall, σ, ϕd, nu = map(1:9) do i
    #     map(outputs) do x
    #         x[i]
    #     end    
    # end
    X, Y, Z, fid, line, zall, σ, ϕd, nu = [[out[i] for out in outputs] for i in 1:9]
end    

function writeasegdat(vonerow::Vector, sfmt::Vector, outfile::String, mode::String)
    #Write single element to file
    open(outfile*".dat", mode) do io
        for (el, fmt) in zip(vonerow, sfmt)
            if isa(el, Array) 
                for iel in el
                print(io, Printf.format(Printf.Format(fmt), iel))    
            end
        else
            print(io, Printf.format(Printf.Format(fmt), el))
            end
        end
    print(io, "\n")
    flush(io)
    close(io)
    end 
end

function writeasegdat(vall::Vector, sfmt::Vector, outfile::String)
    #Write to the file one element at a time across one row
    for i=1:length(vall)
        if i > 1
             mode = "a"
        else 
            mode = "w"
        end
        writeasegdat(vall[i], sfmt, outfile, mode)
    end 
end

function formatasegfield(sfmt::Vector)
    #Create an array with fortran transformned format
    sfmt_fortran = Vector{String}(undef, length(sfmt))
    
    for i =1:length(sfmt)
        temp1 = last(sfmt[i], 1)
        temp1 = uppercase(temp1)
        temp2 = strip(sfmt[i], '%')
    
        fortran_format = string(temp1, temp2)
        fortran_format = chop(fortran_format)
    
        sfmt_fortran[i] = fortran_format
    end
    return sfmt_fortran
end

function writeasegdfnfromonerow(vonerow::Vector, channel_names::Vector, sfmt::Vector, outfile::String)
    sfmt_fortran = formatasegfield(sfmt)
    open(outfile*".dfn", "w") do io
        #println(io, "DEFN   ST=RECD,RT=COMM;RT:A4;COMMENTS:A76") #Geosoft
        for (i, el) in enumerate(vonerow)
            pre = isa(el, Array) ? string(length(el)) : ""
            println(io, "DEFN $i ST=RECD,RT=; $(channel_names[1][i]): "*pre*"$(sfmt_fortran[i]): NULL=-9999999.99, UNITS=$(channel_names[2][i]), LONGNAME=$(channel_names[3][i])")
        end
        println(io, "DEFN $(length(channel_names[1])+1) ST=RECD,RT=; END DEFN")
        flush(io)
        close(io)
    end
end 
 
function writeasegdfn(vall::Vector, channel_names::Vector, sfmt::Vector, outfile::String)
    writeasegdfnfromonerow(vall[1], channel_names::Vector, sfmt::Vector, outfile::String)
end

function readxyzrhoϕ(linenum::Int, nlayers::Int; pathname="")
    # get the rhos
    fnameρ = joinpath(pathname, "rho_avg_line_$(linenum)_summary_xyzrho.txt")
    A = readlargetextmatrix(fnameρ)
    ρavg = reshape(A[:,4], nlayers, :)
    ρlow, ρmid, ρhigh =  map(["low", "mid", "hi",]) do lstring
        fnameρ = joinpath(pathname, "rho_"*lstring*"_line_$(linenum)_summary_xyzrho.txt")
        B = readlargetextmatrix(fnameρ)
        ρ = reshape(B[:,4], nlayers,:)
    end
    # get the X, Y, Z
    X, Y = map(i->A[1:nlayers:end,i], (1:2))
    Zfirstsounding = A[1:nlayers,3] # this is height, does not include surface topo
    zall = getzall(Zfirstsounding)
    Z = A[1:nlayers:end,3] .+ zall[1] # this is topo height
    # get the phid
    ϕmean, ϕsdev =  map(["mean", "sdev"]) do lstring
        fnameϕ = joinpath(pathname, "phid_"*lstring*"_line_$(linenum)_summary.txt")
        ϕ = readlargetextmatrix(fnameϕ)
    end
    X, Y, Z, zall, ρlow, ρmid, ρhigh, ρavg, ϕmean, ϕsdev
end

function readxyzrhoϕnu(linenum::Int, nlayers::Int; pathname="")
    out = readxyzrhoϕ(linenum, nlayers; pathname)
    nulow, numid, nuhigh =  map(["low", "mid", "high",]) do lstring
        fnameρ = joinpath(pathname, "nu_"*lstring*"_line_$(linenum)_summary.txt")
        B = readlargetextmatrix(fnameρ)'
    end
    out..., nulow, numid, nuhigh
end    

function getprobabilisticoutputs(outputs::AbstractArray)
    X, Y, Z, zall, ρlow, ρmid, ρhigh, ρavg, ϕmean, ϕsdev = [[out[i] for out in outputs] for i in 1:10]
end    

function getzall(zheights)
    # first get extendfrac r
    numerator = diff(zheights[2:end])
    denom = diff(zheights[1:end-1])
    r = (denom'*denom)\(denom'*numerator) # overkill but I love least squares
    # now for dz
    numerator = -2diff(zheights)
    denom = map(1:length(zheights)-1) do n
        r^(n-1) + r^n
    end
    dz = (denom'*denom)\(denom'*numerator) # also overkill
    zall,  = setupz(0.0, r; dz, n=length(zheights))
    zall
end

function plotmanygrids(σ, X, Y, Z, zcentre; yl=[], xl=[], cmapσdiff="inferno",
        cmapσ="turbo", vmin=-Inf, vmax=Inf, topowidth=1, fontsize=12, spacefactor=5,
        dr=nothing, dz=nothing, plotbinning=true, δ²=1e-3, regtype=:R1, donn=false, hspace=1,
        figsize=(10,10), smallratio=0.1, preferEright=true, delbin=15., 
        titles=fill("", length(σ)), XYprofiles = nothing, drawprofiles=true, 
        makedifffig=false, vmindiff=Inf, vmaxdiff=Inf, difffigsize=(10,3), lidarfile=nothing,
        leftrightxycorners=[], showbulkstatsfigure=false, regionlinewidth=2, regionlinecolor="k")
    @assert !isnothing(dr) # pass as variable as it is used by other functions too       
    if isnothing(dz)
        mins = [zc[1] for zc in zcentre] # depth to first centre
        dz=2*minimum(mins)
    end    
    nsub = length(σ) + 2 # one invisible subplot
    fig, ax = plt.subplots(nsub, 1, gridspec_kw=Dict("height_ratios" => [ones(nsub-2)..., spacefactor*smallratio, smallratio]),
        figsize=figsize)
    
    binby, binvals = getbinby(X, Y, preferEright)
    flipbycoord!(binby, σ, X, Y, Z)
    rmin, rmax = getrangebinextents(binby)
    r, m, sd = binbycoord(rmin, rmax, delbin, binby, binvals)
    coord_mle = getsmoothline(m, sd; δ², regtype)
    # either of x, y are means, either of xr, yr are the fit
    x, y, xr, yr = get_x_y(r, m, coord_mle, preferEright)
    # now if a lidarfile is provided, get elevation at xr, yr
    if !isnothing(lidarfile)
        zr = getlidarheight(lidarfile, [xr[:]';yr[:]'])
    end
    plotbinning && plotbinningresults(X, Y, x, y, xr, yr)
    outmap = map(zip(σ, X, Y, Z, zcentre)) do (s, xx, yy, topo, zc)
        id = getclosestidx([xx';yy'], xr', yr', showinfo=false)
        if !isnothing(XYprofiles)
            @assert size(XYprofiles, 1) == 2
            idxatXY = getclosestidx([xx';yy'], XYprofiles[1,:]', XYprofiles[2,:]', showinfo=false)
            zatXY =  zc
            satXY = s[:,idxatXY]
        else
            zatXY = nothing 
            satXY = nothing
        end
        topouse = isnothing(lidarfile) ? topo[id] : zr
        img, gridr, gridz, topofine, R = makegrid(s[:,id], xr, yr, topouse; donn, dr, zall=zc, dz)
        img, gridr, gridz, topofine, R, zatXY, satXY
    end
    img, gridr, gridz, topofine, R, zatXY, satXY = [[out[i] for out in outmap] for i in 1:7]
    if (isinf(vmin) || isinf(vmax))
        vmin, vmax = extrema(reduce(vcat, [[extrema(s)...] for s in σ]))
    end
    if !isnothing(XYprofiles)
        @assert size(XYprofiles, 1) == 2
        idx_ = getclosestidx([xr'; yr'], XYprofiles[1,:]', XYprofiles[2,:]', showinfo=false)
        Rtoplotat = [cumulativelinedist(xr, yr)[id_] for id_ in idx_]
    end
    bulkstat_mean = !isempty(leftrightxycorners) ? zeros(length(σ)) : zeros(0)
    imhandle = map(zip(ax, img, gridr, gridz, topofine, titles, 1:length(σ))) do (
        ax_, img_, gridr_, gridz_, topofine_, ti, i) 
        imhandle_ = ax_.imshow(img_, extent=[gridr_[1], gridr_[end], gridz_[end], gridz_[1]]; 
            cmap=cmapσ, aspect="auto", vmin, vmax)
        ax_.plot(gridr_, topofine_, linewidth=topowidth, "-k")
        if (!isnothing(XYprofiles) && drawprofiles)
            idx, = nn(KDTree(gridr_'), [Rtoplotat...]')
            plotprofile(ax_, idx, topofine_, gridr_)
        end    
        isempty(ti) || ax_.set_title(ti)
        if !isempty(leftrightxycorners)
            bulkstat_mean[i] = getbulkstatfromimageregion(img_, gridr_, gridz_; 
                leftrightxycorners, showfig=showbulkstatsfigure)
        end
        imhandle_
    end
    map(1:nsub-3) do i
        ax[i].sharex(ax[i+1])
        ax[i].sharey(ax[i+1])
    end
    [a.tick_params(labelbottom=false) for a in ax[1:end-3]]
    [a.set_ylabel("Height m") for a in ax[1:end-2]]
    ax[end-2].set_xlabel("Distance m")
    ax[end-1].axis("off")
    !isempty(xl) && ax[1].set_xlim(xl)
    isempty(yl) && (yl = extrema(reduce(vcat, gridz)))
    ax[1].set_ylim(yl)
    cb = fig.colorbar(imhandle[end], cax=ax[end], orientation="horizontal")
    cb.set_label("Log₁₀ S/m", labelpad=0)
    nicenup(fig, fsize=fontsize, h_pad=0)
    fig.subplots_adjust(hspace = all(isempty.(titles)) ? 0 : hspace)

    if makedifffig
        difffig(img, gridr[1], gridz[1], topofine[1], difffigsize, vmindiff, vmaxdiff, cmapσdiff, topowidth, yl, 
            fontsize, leftrightxycorners, regionlinecolor, regionlinewidth)
    end    

    xr, yr, ax, zatXY, satXY, bulkstat_mean # return easting northing of grid and figure axes
end

function getbulkstatfromimageregion(imgin, gridr, gridz; leftrightxycorners=[], central_tendency=:mean, showfig=false)
    @assert !isempty(leftrightxycorners)
    x1, elev1, x2, elev2 = leftrightxycorners
    @assert x1<x2
    @assert elev1>elev2
    idxr = x1 .< gridr .< x2
    idxz = elev1 .> gridz .> elev2
    # won't work for edge pixel at right or bottom!!
    img = copy(imgin)
    mpart = vec(img[idxz[1:end-1], idxr[1:end-1]])
    # @info sum(.!isnan.(mpart))
    m = eval(central_tendency)(mpart[.!isnan.(mpart)])
    img[idxz[1:end-1], idxr[1:end-1]] .= m
    if showfig
        f, ax = plt.subplots(1,1)
        ax.imshow(img, extent=[gridr[1], gridr[end], gridz[end], gridz[1]];aspect="auto")
    end
    m
end

function difffig(img, gridr, gridz, topofine, difffigsize, vmindiff, vmaxdiff, cmapσdiff, topowidth, 
        yl, fontsize, leftrightxycorners, color, linewidth)
    fig, ax = plt.subplots(1,1, figsize=difffigsize)
    sdimage = std(img)
    if (isinf(vmindiff) || isinf(vmaxdiff))
        vmindiff, vmaxdiff = extrema(.!isnan.(sdimage))
    end
    imhandle = ax.imshow(std(img), extent=[gridr[1], gridr[end], gridz[end], gridz[1]]; 
            cmap=cmapσdiff, aspect="auto", vmin=vmindiff, vmax=vmaxdiff)
    ax.plot(gridr, topofine, linewidth=topowidth, "-k")
    if !isempty(leftrightxycorners)
        x1, elev1, x2, elev2 = leftrightxycorners
        ax.add_patch(matplotlib.patches.Rectangle((x1,elev1), x2-x1, (elev2-elev1); edgecolor=color, fill=false, linewidth))
    end
    ax.set_ylim(yl)
    cb = fig.colorbar(imhandle)
    cb.set_label("Log₁₀ S/m")
    ax.set_ylabel("Height m")
    ax.set_xlabel("Distance m")
    nicenup(fig, fsize=fontsize)
end

function getbinby(X, Y, preferEright)
    binby, binvals = X, Y
    if !preferEright
        binby, binvals = Y, X
    end
    binby, binvals
end

function plotbinningresults(X, Y, x, y, xr, yr, plotmean=false)
    _, ax = plt.subplots(1, 1)
    for i in eachindex(X)
        ax.plot(X[i],Y[i],"-", color="gray") #label=split(fnames[i],"/")[2])
    end
    plotmean && ax.plot(x, y,"-k", label="mean path")
    ax.plot(xr, yr, linewidth=1, color="r", label="mle path")
    ax.legend()
    ax.grid()
    ax.set_aspect(1)
end    

function get_x_y(r, m, coord_mle, preferEright)
    if preferEright
        x, y = r, m
        xr, yr = x, coord_mle 
    else
        x, y = m, r
        xr, yr = coord_mle, y
    end  
    x, y, xr, yr
end    

function flipbycoord!(coordsarray, stufftoflip...) # slurp
    for (i, coords) in enumerate(coordsarray)
        if coords[end]<coords[1]
            _ = map(stufftoflip) do x 
                if size(x[i], 2) == 1
                    reverse!(x[i])    
                else
                    reverse!(x[i], dims=2)
                end
            end    
        end
    end
end

function getrangebinextents(XX)
    rmin = maximum([x[1] for x in XX])
    rmax = minimum([x[end] for x in XX])
    rmin, rmax
end

function binbycoord(rmin, rmax, delbin, binby, binvals,)
    r = range(rmin, rmax, step=delbin)
    m, sd  = map(x->zeros(length(r)-1), 1:2)
    for i in 2:length(r)
        s, s2, n = 0., 0., 0
        for (bin, val) in zip(binby, binvals)
            # below ensures closed on both ends
            if i==2
                idx = r[i-1] .<= bin .<= r[i]
            else
                idx = r[i-1] .< bin .<= r[i]
            end
            if !isempty(idx)
                s  += sum(val[idx])
                s2 += sum(val[idx].^2)
                n  += sum(idx)
            end    
        end
        @assert n>2 "nothing to bin, make delbin larger than $delbin"
        m[i-1]  = s/n
        varest = s2/n - m[i-1]^2
        varest < 0.25 && (varest = 1) # variance estimates can be wonky if coincident points, few points etc.
        sd[i-1] = sqrt(varest)
    end
    (r[1:end-1]+r[2:end])/2, m, sd        
end

function readwell(fname, skipstart; lidarfile=nothing)
    # in Ross' .con format
    # skip some lines and then must have format in the skipped lines
    # bore: thisborename
    # and then
    # depth mS/m
    # returns converted to log10 S/m
    io = open(fname)
    name = ""
    X, Y, Z = 0., 0., 0.
    for (i, str) in enumerate(eachline(io))
        s = split(str, ":")
        lowercase(s[1]) == "bore" && (name = s[end])
        lowercase(s[1]) == "easting"   && (X = parse(Float64, s[end]))
        lowercase(s[1]) == "northing"  && (Y = parse(Float64, s[end]))
        if lowercase(s[1]) == "elevation" 
            if isnothing(lidarfile)
                Z = parse(Float64, s[end])
            else # get it from a lidar point cloud
                Z = getlidarheight(lidarfile,[X;Y])
            end    
        end            
        i==skipstart && break
    end    
    @info name
    zc_rho = readdlm(fname; skipstart)
    zc_rho[:,2] = 3 .-log10.(zc_rho[:,2]) # log10 ohm m
    name, X, Y, Z, zc_rho
end

function findMAE(soundings::Vector{T}, opt_in::OptionsStat, Xwell, Ywell, zc_rho; burninfrac=0.5) where T<:Sounding
    opt = deepcopy(opt_in)
    zcentre = opt.xall
    Mwell = blocktomodel(zcentre, zc_rho)
    m = getclosestensemble(soundings, opt, Xwell, Ywell; burninfrac)
    findMAE(m, vec(Mwell))
end

function blocktomodel(zcentre, zc_rho)
    zboundaries = zcentertoboundary(zcentre)
    Mwell = block1Dvalues([zc_rho[:,2]], zc_rho[:,1], [zboundaries[1:end-1] zboundaries[2:end]], :median)[:]
    Mwell = vec([Mwell;Mwell[end,:]']) # dummy last cell in depth
end

function getclosestensemble(soundings::Vector{T}, opt::OptionsStat, Xwell, Ywell; burninfrac=0.5) where T<:Sounding
    idx = getclosestidx(Xwell, Ywell, soundings)
    opt.fdataname = soundings[idx].sounding_string*"_"
    m = assembleTat1(opt, :fstar; burninfrac, temperaturenum=1)
end

function findMAE(mm::Vector{T}, Mwell::Vector{S}) where T<:Array where S<:Real
    absdevs = [abs.(vec(m) - Mwell) for m in mm]
    nanmean(reduce(hcat, absdevs), 2) # mean of ndepths × namples in nsamples dirn
end

function makeblockedwellimage(wellarray, zall, xr, yr; distblank=50, dr=nothing, donn=false, displaydistanceaway=true, lidarfile=nothing )
    # xr, yr are the line path along which to find closest well index
    @assert !isnothing(dr)
    @assert mod(distblank, dr) == 0
    wellname, Xwell, Ywell, Zwell, z_rho_well = [[well[i] for well in wellarray] for i in 1:5]
    zboundaries = zcentertoboundary(zall)
    Mwell = reduce(hcat, map(z_rho_well) do zρ
                block1Dvalues([zρ[:,2]], zρ[:,1], [zboundaries[1:end-1] zboundaries[2:end]], :mean)'
    end)
    Mwell = [Mwell;Mwell[end,:]'] # dummy last cell in depth
    # idx is length of xr, yr
    idx, _ = getclosestidxanddist([Xwell';Ywell'], xr', yr')
    if !isnothing(lidarfile)
        zr = getlidarheight(lidarfile, [xr[:]';yr[:]'])
    end
    topouse = isnothing(lidarfile) ? Zwell[idx] : zr
    Mclosest = Mwell[:,idx] # this needs to be plotted on image of line with coordinates xr, yr
    # this is inline distance along line
    rr = cumulativelinedist(xr, yr)
    # interpolate linearly as usual onto line with xr, yr coordinates with depth and line distance
    img, gridr, gridz, topofine = makegrid(Mclosest, xr, yr, topouse; donn, dr, zall, dz=zall[1]*2)
    # idxclosest is length of wellarray
    idxclosest, distanceaway = getclosestidxanddist([xr'; yr'], Xwell', Ywell')
    displaydistanceaway && [@printf("WELL: %s DISTANCE: %.2f m\n", w, d) for (w,d) in zip(wellname, distanceaway)]
    # NaN out farther than distblank
    gridcenter = 0.5(gridr[1:end-1]+gridr[2:end])
    dist = reduce(hcat, [abs.(gridcenter .- rr[i]) for i in idxclosest])
    distflag = vec(reduce(&, dist .> distblank, dims=2))
    img[:,distflag] .= NaN
    # now find changes in the flag to mark well location and refine idxclosest to conform to grid
    c = findall(distflag[:] .== 0)
    idxclosest = zeros(Int, length(wellarray))
    idxclosest[1] = c[1]
    count = 1
    for i in 2:length(c)
        if c[i-1] != c[i]-1
            count += 1
            idxclosest[count] = c[i]
        end
    end
    xcoords, elevcoords, depth = rectangulatewells(gridcenter, zall, Mwell, distblank, topofine, idxclosest, dr)
    img, gridr, gridz, xcoords, elevcoords, depth
end    

function rectangulatewells(gridcenter, zall, Mwell, distblank, topouse, idxclosest, dr)
    xcoords = gridcenter[idxclosest] .- dr/2
    idx_z_start = [findfirst(.!isnan.(m)) for m in eachcol(Mwell)] .-1
    elevcoords = topouse[idxclosest] - zall[idx_z_start] 
    idx_z_end = [findlast(.!isnan.(m)) for m in eachcol(Mwell)] 
    depths = map(zip(idx_z_start, idx_z_end)) do (izstart, izend)
        depth = zall[izend] - zall[izstart]
    end
    xcoords, elevcoords, depths
end

function plotwelloutline(ax, img, xcoords, elevcoords, depths, gridr, gridz, vmin, vmax, distblank; cmap="turbo", color="k", linewidth=0.5)
    # plot the well outline on axis
    ax.imshow(-img; extent=[gridr[1], gridr[end], gridz[end], gridz[1]], cmap, vmin, vmax)
    map(zip(xcoords, elevcoords, depths)) do (x, y, Δy)
        ax.add_patch(matplotlib.patches.Rectangle((x,y), 2distblank, -Δy; edgecolor=color, fill=false, linewidth))
    end 
    ax.set_aspect("auto")
end

function plotwelloutline(ax::Array, img, xcoords, elevcoords, depths, gridr, gridz, vmin, vmax, distblank; cmap="turbo", color="k", linewidth=0.5)
    # plot into each axis, the outline
    for a in ax
        plotwelloutline(a, img, xcoords, elevcoords, depths, gridr, gridz, vmin, vmax, distblank; cmap, color, linewidth)
    end
end

function getlidarheight(lidarheightfile::String, xy)
    A = readlargetextmatrix(lidarheightfile)
    # read X, Y, mAHDevery 10 m
    kdtree = KDTree(A[:,1:2]')
    if size(xy, 2) == 2     
        idxs, = nn(kdtree,xy')
    else
        idxs, = nn(kdtree,xy)
    end
    A[idxs,3]
end    

function plotblockedwellonimages(ax, wellarray, zall, xr, yr; donn=false,
        vmin=nothing, vmax=nothing, dr=15, distblank=4dr, cmap="turbo", color="k", linewidth=0.5, lidarfile=nothing)
        @assert !isnothing(vmin)
        @assert !isnothing(vmax)
    img, gridr, gridz, xcoords, elevcoords, depths = makeblockedwellimage(wellarray, zall, xr, yr; distblank, dr, 
        lidarfile)
    plotwelloutline(ax, img, xcoords, elevcoords, depths, gridr, gridz, vmin, vmax, distblank; cmap, color, linewidth)
    wnames = [String(w[1]) for w in wellarray]
    annotatewells(ax[end], wnames, xcoords, elevcoords, depths, distblank)
    
end

function annotatewells(ax, wnames, xcoords, elevcoords, depths, distblank)
    map(zip(xcoords, elevcoords, depths, wnames)) do (x, y, Δy, wn)
        ax.text(x +2distblank, y-Δy, wn, fontsize=8)
    end
end

end # module CommonToAll
