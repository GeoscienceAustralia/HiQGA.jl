using TransD_GP, PyPlot, StatsBase, Statistics, LinearAlgebra

mutable struct Line<:Operator
    n::Nothing
end

function plot_posterior(L::Line,
                        optns::TransD_GP.OptionsNonstat,
                        opt::TransD_GP.OptionsStat;
    nbins = 50,
    burninfrac=0.5,
    isns=true,
    qp1=0.05,
    qp2=0.95,
    cmappdf = "RdYlBu_r",
    figsize=(10,5),
    pdfnormalize=false)
    himage_ns, edges_ns, CI_ns = make1Dhist(optns, burninfrac=burninfrac, nbins = nbins, qp1=qp1, qp2=qp2, pdfnormalize=pdfnormalize)
    himage, edges, CI = make1Dhist(opt, burninfrac=burninfrac, nbins = nbins, qp1=qp1, qp2=qp2, pdfnormalize=pdfnormalize)
    f,ax = plt.subplots(1,2, sharey=true, figsize=figsize)
    xall = opt.xall
    im1 = ax[1].imshow(himage_ns, extent=[edges_ns[1],edges_ns[end],xall[end],xall[1]], aspect="auto", cmap=cmappdf)
    cb1 = colorbar(im1, ax=ax[1])
    cb1.ax.set_xlabel("pdf ns")
    ax[1].grid()
    im2 = ax[2].imshow(himage, extent=[edges[1],edges[end],xall[end],xall[1]], aspect="auto", cmap=cmappdf)
    ax[2].grid()
    cb2 = colorbar(im2, ax=ax[2])
    cb2.ax.set_xlabel("pdf stationary")
    ax[1].plot(CI_ns, xall[:], linewidth=2, color="g")
    ax[1].plot(CI_ns, xall[:], linewidth=2, color="k", linestyle="--")
    ax[2].plot(CI, xall[:], linewidth=2, color="g")
    ax[2].plot(CI, xall[:], linewidth=2, color="k", linestyle="--")
    ax[1].set_xlabel(L"\log_{10} \rho")
    ax[1].set_ylabel("depth (m)")
    ax[2].set_xlabel(L"\log_{10} λ")
end

function make1Dhist(opt::TransD_GP.Options;
                burninfrac = 0.5,
                nbins = 50,
                rhomin=Inf,
                rhomax=-Inf,
                qp1=0.05,
                qp2=0.95,
                pdfnormalize=false)
    M = assembleTat1(opt, burninfrac=burninfrac)
    if (rhomin == Inf) && (rhomax == -Inf)
        for (i,mm) in enumerate(M)
            rhomin_mm = minimum(mm)
            rhomax_mm = maximum(mm)
            rhomin_mm < rhomin && (rhomin = rhomin_mm)
            rhomax_mm > rhomax && (rhomax = rhomax_mm)

        end
        if typeof(opt) == TransD_GP.OptionsStat
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
        if typeof(opt) == TransD_GP.OptionsStat
            himage[ilayer,:] = fit(Histogram, [0.5log10.(m[ilayer]) for m in M], edges).weights
        else
            himage[ilayer,:] = fit(Histogram, [m[ilayer] for m in M], edges).weights
        end
        himage[ilayer,:] = himage[ilayer,:]/sum(himage[ilayer,:])/(diff(edges)[1])
        pdfnormalize && (himage[ilayer,:] = himage[ilayer,:]/maximum(himage[ilayer,:]))
        if typeof(opt) == TransD_GP.OptionsStat
            CI[ilayer,:] = [quantile([0.5log10.(m[ilayer]) for m in M],(qp1, qp2))...]
        else
            CI[ilayer,:] = [quantile([m[ilayer] for m in M],(qp1, qp2))...]
        end
    end
    return himage, edges, CI
end

function get_posterior(opt_in::TransD_GP.Options, stat::Symbol; decimate=1)
    TransD_GP.history(opt_in, stat=stat)[1:decimate:end]
end

function assembleTat1(optin::TransD_GP.Options, stat::Symbol; burninfrac=0.5)
    @assert (stat == :fstar || stat == :x_ftrain)
    isns = checkns(optin)
    @assert 0.0<=burninfrac<=1.0
    Tacrosschains = gettargtemps(optin)
    iters = TransD_GP.history(optin, stat=:iter)
    start = round(Int, length(iters)*burninfrac)
    start == 0 && (start = 1)
    @info "obtaining models $(iters[start]) to $(iters[end])"
    nmodels = sum((Tacrosschains[start:end,:] .== 1))
    mat1 = Array{Array{Float64}, 1}(undef, nmodels)
    opt = deepcopy(optin)
    imodel = 0
    for ichain in 1:size(Tacrosschains, 2)
        at1idx = findall(Tacrosschains[:,ichain].==1).>= start
        ninchain = sum(at1idx)
        @info "chain $ichain has $ninchain models"
        ninchain == 0 && continue
        opt.fstar_filename = "models_"*isns*opt.fdataname*"_$ichain.bin"
        opt.x_ftrain_filename = "points_"*isns*opt.fdataname*"_$ichain.bin"
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

get_misfit(m::TransD_GP.ModelNonstat, opt::TransD_GP.Options, L::Line) = 0.0

export Line
