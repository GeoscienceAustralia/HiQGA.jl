using TransD_GP, PyPlot, StatsBase, Statistics, LinearAlgebra

mutable struct Line<:Operator
    n::Nothing
end

function plot_posterior(L::Line, opt_in::TransD_GP.Options;
    nbins = 50,
    burninfrac=0.5,
    isns=true,
    qp1=0.05,
    qp2=0.95,
    cmapcdf = "RdYlBu_r",
    cmappdf = "inferno_r",
    figsize=(10,5),
    pdfnormalize=false)
    M = assembleTat1(opt_in, burninfrac=burninfrac, isns=isns)
    xall = opt_in.xall
    f,ax = plt.subplots(1,2, sharex=true, sharey=true, figsize=figsize)
    rhomin, rhomax = Inf, -Inf
    for (i,mm) in enumerate(M)
        rhomin_mm = minimum(mm)
        rhomax_mm = maximum(mm)
        rhomin_mm < rhomin && (rhomin = rhomin_mm)
        rhomax_mm > rhomax && (rhomax = rhomax_mm)
    end
    edges = LinRange(rhomin, rhomax, nbins+1)
    himage = zeros(Float64, length(M[1]), nbins)
    cumhimage = zeros(Float64, length(M[1]), nbins)
    CI = zeros(Float64, length(M[1]), 2)
    for ilayer=1:size(xall,2)
        himage[ilayer,:] = fit(Histogram, [m[ilayer] for m in M], edges).weights
        himage[ilayer,:] = himage[ilayer,:]/sum(himage[ilayer,:])/(diff(edges)[1])
        cumhimage[ilayer,:] = cumsum(himage[ilayer,:])/sum(himage[ilayer,:])
        pdfnormalize && (himage[ilayer,:] = himage[ilayer,:]/maximum(himage[ilayer,:]))
        CI[ilayer,:] = [quantile([m[ilayer] for m in M],(qp1, qp2))...]
    end
    im1 = ax[1].imshow(himage, extent=[edges[1],edges[end],xall[end],xall[1]], aspect="auto",
            cmap=cmappdf)
    cb1 = colorbar(im1, ax=ax[1])
    cb1.ax.set_xlabel("pdf")
    ax[1].grid()
    im2 = ax[2].imshow(cumhimage, extent=[edges[1],edges[end],xall[end],xall[1]], aspect="auto",
            cmap=cmapcdf)
    ax[2].grid()
    cb2 = colorbar(im2, ax=ax[2])
    cb2.ax.set_xlabel("cdf")
    ax[1].plot(CI, xall[:], linewidth=2, color="g")
    ax[1].plot(CI, xall[:], linewidth=2, color="k", linestyle="--")
    ax[2].plot(CI, xall[:], linewidth=2, color="g")
    ax[2].plot(CI, xall[:], linewidth=2, color="k", linestyle="--")
    ax[1].set_xlabel(L"\log_{10} \rho")
    ax[1].set_ylabel("depth (m)")
    ax[2].set_xlabel(L"\log_{10} \rho")
end

function get_posterior(opt_in::TransD_GP.Options, stat::Symbol; decimate=1)
    TransD_GP.history(opt_in, stat=:fstar)[1:decimate:end]
end

function assembleTat1(optin::TransD_GP.Options; burninfrac=0.5, isns=true)
    @info "ns is $isns"
    ns = "ns"
    isns || (ns="")
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
        @info "chain $ichain"
        at1idx = findall(Tacrosschains[:,ichain].==1).>= start
        ninchain = sum(at1idx)
        ninchain == 0 && continue
        opt.fstar_filename = "models_"*ns*opt.fdataname*"_$ichain.bin"
        mat1[imodel+1:imodel+ninchain] .= TransD_GP.history(opt, stat=:fstar)[at1idx]
        imodel += ninchain
    end
    mat1
end

function gettargtemps(opt_in::TransD_GP.Options)
    nchains = length(filter( x -> occursin(r"misfits_ns.*bin", x), readdir(pwd()) )) # my terrible regex
    @info "Number of chains is $nchains"
    # now look at any chain to get how many iterations
    opt = deepcopy(opt_in)
    costs_filename = "misfits_ns"*opt.fdataname
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
