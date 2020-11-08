module PlotEM1D
using MCMC_Driver, TransD_GP, PyPlot, JLD, StatsBase, Statistics, CSEM1Dkernels, LinearAlgebra

function get_posterior(opt_in::TransD_GP.Options, opt_EM::MCMC_Driver.EMoptions;
    fdataname="",
    burnin=300000,
    nbins=50,
    rho=[0.],
    alpha=0.025,
    alphacolor = "blue",
    alphalwidth = 0.3,
    qp1=0.01,
    qp2=0.99,
    figsize = (10,6),
    fontsize = 16,
    convfigsize = (8,6),
    kfigsize = (9,6),
    convfontsize = 14,
    cmapcdf = "RdYlBu_r",
    cmappdf = "inferno_r",
    showtitle=false,
    titlestr="",
    topadjust=0.9,
    fixviewextent=[0. 0. 0. 0.],
    rhobgcolor="yellow",
    rhofgcolor="k",
    nlldelta = 100,
    nxbins = 20,
    nftbins = 20,
    nxftfontsize = 16,
    pdfnormalize=false)
    @assert length(fixviewextent) == 4

    M, n, x_ft, iter, misfit, T0loc = get_posterior(opt_in, fdataname=fdataname)

    from = findfirst(iter.>burnin)
    f1 = figure(figsize=convfigsize)
    s1 = subplot(311)
    plot(iter,n)
    s1[:set_ylim](opt_in.nmin, opt_in.nmax)
    ylabel("# train")
    s2 = subplot(312, sharex=s1)
    plot(iter, misfit)
    ylabel("-log L")
    s3 = subplot(313, sharex=s1)
    plot(iter, T0loc)
    s3[:set_yticks](collect(minimum(T0loc):2:maximum(T0loc)))
    s2[:set_ylim](mean(misfit[from:end])-nlldelta, mean(misfit[from:end])+nlldelta)
    s3[:set_xlim](burnin+1, iter[end])
    ylabel("T=1 location")
    xlabel("iterations")
    MCMC_Driver.nicenup(gcf(), fsize=convfontsize)
    if showtitle
        titlestr == "" && (titlestr = fdataname)
        f1[:suptitle](titlestr, fontsize=convfontsize)
        f1[:subplots_adjust](top=topadjust)
    end

    xall = opt_in.xall
    f,ax = plt[:subplots](1,2, sharex=true, sharey=true, figsize=figsize)
    rhomin, rhomax = 0.0, 0.0
    for (i,mm) in enumerate(M)
        if i >= from
            # ax[1][:step](mm, xall[:], linewidth=alphalwidth, color=alphacolor, alpha=alpha)
            #plot(x_ft[i][1:n[i],2], x_ft[i][1:n[i],1], ".", markersize=20)
            # title("model $(iter[i])")
        end
        rhomin_mm = min(mm...)
        rhomax_mm = max(mm...)
        rhomin_mm < rhomin && (rhomin = rhomin_mm)
        rhomax_mm > rhomax && (rhomax = rhomax_mm)
    end

    # if rho != [0.]
    #     ax[1][:step](log10.(rho[2:end]), opt_EM.z[2:end], linewidth=3, color=rhobgcolor)
    #     ax[1][:step](log10.(rho[2:end]), opt_EM.z[2:end], linewidth=2, color=rhofgcolor)
    # end
    # ax[1][:grid]()
    # ax[1][:invert_yaxis]()

    edges = LinRange(rhomin, rhomax, nbins+1)
    himage = zeros(Float64, length(M[1]), nbins)
    cumhimage = zeros(Float64, length(M[1]), nbins)
    CI = zeros(Float64, length(M[1]), 2)
    z = opt_EM.z
    for ilayer=1:length(M[1])
        himage[ilayer,:] = (fit(Histogram, [m[ilayer] for m in M[from:end]], edges).weights)
        himage[ilayer,:] = himage[ilayer,:]/sum(himage[ilayer,:])/(diff(edges)[1])
        cumhimage[ilayer,:] = cumsum(himage[ilayer,:])/sum(himage[ilayer,:])
        pdfnormalize && (himage[ilayer,:] = himage[ilayer,:]/maximum(himage[ilayer,:]))
        CI[ilayer,:] = [quantile([m[ilayer] for m in M[from:end]],(qp1, qp2))...]
    end
    im1 = ax[1][:imshow](himage, extent=[edges[1],edges[end],xall[end],xall[1]], aspect="auto",
            cmap=cmappdf)
    cb1 = colorbar(im1, ax=ax[1])
    cb1[:ax][:set_xlabel]("pdf")
    ax[1][:grid]()
    im2 = ax[2][:imshow](cumhimage, extent=[edges[1],edges[end],xall[end],xall[1]], aspect="auto",
            cmap=cmapcdf)
    ax[2][:grid]()
    cb2 = colorbar(im2, ax=ax[2])
    cb2[:ax][:set_xlabel]("cdf")
    if rho != [0.]
        ax[1][:step](log10.(rho[2:end]), opt_EM.z[2:end], linewidth=3, color=rhobgcolor)
        ax[1][:step](log10.(rho[2:end]), opt_EM.z[2:end], linewidth=2, color=rhofgcolor)
        ax[2][:step](log10.(rho[2:end]), opt_EM.z[2:end], linewidth=3, color=rhobgcolor)
        ax[2][:step](log10.(rho[2:end]), opt_EM.z[2:end], linewidth=2, color=rhofgcolor)
    end
    ax[1][:plot](CI, xall[:], linewidth=2, color="g")
    ax[1][:plot](CI, xall[:], linewidth=2, color="k", linestyle="--")
    ax[2][:plot](CI, xall[:], linewidth=2, color="g")
    ax[2][:plot](CI, xall[:], linewidth=2, color="k", linestyle="--")
    # ax[2][:plot](CI, xall[:], "--k")
    ylim(opt_EM.zRx, xall[end])
    gca()[:invert_yaxis]()
    ax[1][:set_xlabel](L"\log_{10} \rho")
    ax[1][:set_ylabel]("depth (m)")
    ax[2][:set_xlabel](L"\log_{10} \rho")
    if !(fixviewextent == [0. 0. 0. 0.])
        ax[1][:set_xlim](fixviewextent[1:2])
        ax[1][:set_ylim](fixviewextent[3:4])
        ax[1][:invert_yaxis]()
    end
    MCMC_Driver.nicenup(f, fsize=fontsize)
    if showtitle
        titlestr == "" && (titlestr = fdataname)
        f[:suptitle](titlestr, fontsize=fontsize)
        f[:subplots_adjust](top=topadjust)
    end

    x, ft = zeros(sum(n[from:end])), zeros(sum(n[from:end]))
    nlast = 0
    for (i, xt) in enumerate(x_ft[from:end])
        x[nlast+1:nlast+n[from-1+i]]  = xt[1:n[from-1+i],1]
        ft[nlast+1:nlast+n[from-1+i]] = xt[1:n[from-1+i],2]
        nlast = nlast+n[from-1+i]
    end
    f2, ax2 = plt[:subplots](2,2, figsize=kfigsize)
    edgesx = LinRange(opt_in.xbounds[1], opt_in.xbounds[2], nxbins+1)
    edgesrho = LinRange(opt_in.fbounds[1], opt_in.fbounds[2], nftbins+1)
    h = fit(Histogram, (x, ft), (edgesx, edgesrho)).weights
    im3 = ax2[1][:imshow](h, extent=[edgesrho[1],edgesrho[end],edgesx[end],edgesx[1]], aspect="auto", vmin=0.0, vmax=2*max(h...))
    ax2[1][:set_ylabel]("depth m")
    ax2[1][:set_xlabel](L"\log_{10}\rho")
    rhist = sum(h,dims=1)[:]
    rhist = rhist/sum(rhist)/diff(edgesrho)[1]
    ax2[2][:bar](0.5*(edgesrho[2:end]+edgesrho[1:end-1]), rhist, width=diff(edgesrho)[1])
    ax2[2][:plot]([edgesrho[1], edgesrho[end]], 1/(opt_in.fbounds[end]-opt_in.fbounds[1])*[1,1],"--k")
    ax2[2][:set_ylim](0, 3*maximum(rhist))
    ax2[2][:set_xlabel](L"\log_{10}\rho")
    ax2[2][:set_ylabel]("pdf")
    xhist = sum(h,dims=2)[:]
    xhist = xhist/sum(xhist)/diff(edgesx)[1]
    ax2[3][:barh](0.5*(edgesx[2:end]+edgesx[1:end-1]), xhist, height=diff(edgesx)[1])
    ax2[3][:plot](1/(opt_in.xbounds[end]-opt_in.xbounds[1])*[1,1], [edgesx[1], edgesx[end]],"--k")
    ax2[3][:set_xlim](0, 3*maximum(xhist))
    ax2[3][:set_ylim](ax2[1][:get_ylim]())
    ax2[3][:set_ylabel]("depth m")
    ax2[3][:set_xlabel]("pdf")
    # ax2[3][:xaxis][:set_label_position]("top")
    ax2[2][:get_shared_y_axes]()[:join](ax2[1], ax2[3])
    ax2[2][:get_shared_x_axes]()[:join](ax2[1], ax2[2])
    ax2[2][:set_xlim](ax2[1][:get_xlim]())
    k = fit(Histogram, n[from:end], (opt_in.nmin-0.5:opt_in.nmax+0.5)).weights
    k = k/sum(k)
    @show mean(k)
    ax2[4][:bar](opt_in.nmin:opt_in.nmax, k, width=1)
    ax2[4][:plot]([opt_in.nmin, opt_in.nmax],1/(opt_in.nmax-opt_in.nmin+1)*[1,1],"--k")
    ax2[4][:set_xlabel]("# training points")
    ax2[4][:set_ylim](0, 3*max(k...))
    ax2[4][:set_ylabel]("pdf")
    ax2[4][:set_xlim](opt_in.nmin, opt_in.nmax)
    MCMC_Driver.nicenup(gcf(), fsize=nxftfontsize)
    savefig(fdataname*"_k.png", dpi=300)

    figure(f1[:number])
    savefig(fdataname*"_conv.png", dpi=300)
    figure(f[:number])
    savefig(fdataname*"_post.png", dpi=300)
    return ax[1], im1, cb1, ax[2], im2, cb2
end

function get_posterior(opt_in::TransD_GP.Options; fdataname="", decimate=1)

    M = TransD_GP.history(opt_in, stat=:fstar)
    n = TransD_GP.history(opt_in, stat=:nodes)
    x_ft = TransD_GP.history(opt_in, stat=:x_ftrain)
    iter = TransD_GP.history(opt_in, stat=:iter)
    misfit = load("misfit_T0_"*fdataname*".jld", "misfit")
    T0loc = load("misfit_T0_"*fdataname*".jld", "T0loc")

    return M[1:decimate:end], n[1:decimate:end],
    x_ft[1:decimate:end], iter[1:decimate:end], misfit[1:decimate:end], T0loc[1:decimate:end]
end

function get_field(mfstar::AbstractArray, d::AbstractArray,
    opt::TransD_GP.Options, opt_EM::MCMC_Driver.EMoptions;
    plotstackerrors=true, plotMLerrors=true, sfactor=1.0)
    @assert length(opt_EM.rho) == length(mfstar) + opt_EM.nfixed

    opt_EM.rho[opt_EM.nfixed+1:end] = 10.0.^mfstar

    Er = CSEM1Dkernels.getCSEM1DanisoHED(opt_EM.freqs,
    opt_EM.rRx, opt_EM.zRx, opt_EM.zTx, opt_EM.z, [opt_EM.rho opt_EM.rho],
    0, RxAzim = opt_EM.RxAzim, TxDip = opt_EM.TxDip)

    σ = zeros(size(Er))
    if plotMLerrors || plotstackerrors
        fig, ax = plt[:subplots](1,2, sharex=true, figsize=(12,8))
    end
    for lfreq = 1:length(opt_EM.freqs)
        nr = sum(.!isnan.(d[:,lfreq]))
        select = .!isnan.(d[:,lfreq])
        r = Er[select,lfreq] - d[select,lfreq]
        σ[:,lfreq] = sqrt(0.5/nr)*norm(r./abs.(d[select,lfreq]), 2)*abs.(d[:,lfreq])
        if plotstackerrors
            ax[1][:errorbar](opt_EM.rRx[select], abs.(d[select,lfreq]), yerr = 2*sfactor*opt_EM.sd[select,lfreq], linestyle="none", marker="x",
            elinewidth=0, capsize=3)
            ax[2][:errorbar](opt_EM.rRx[select], unwrap(angle.(d[select,lfreq])), yerr = 2*sfactor*opt_EM.sd[select,lfreq]./abs.(d[select,lfreq]),
            linestyle="none", marker="x", elinewidth=0, capsize=3)
        end
        if plotMLerrors
            ax[1][:errorbar](opt_EM.rRx[select], abs.(Er[select,lfreq]), yerr = 2*σ[select,lfreq],
            elinewidth=2, color="k")
            ax[2][:errorbar](opt_EM.rRx[select], unwrap(angle.(Er[select,lfreq])), yerr = 2*σ[select,lfreq]./abs.(d[select,lfreq]),
            elinewidth=2, color="k")
        end
    end
    ax[1][:set_yscale]("log")
    ax[1][:grid]()
    ax[2][:grid]()
    return Er, σ
end

function plot_post_fields(d::AbstractArray,
    opt::TransD_GP.Options,
    opt_EM::MCMC_Driver.EMoptions,
    fdataname::String;
    sfactor=1.0, decimate=1000, burnin=0, alpha=0.1, figsize=(10,7), fig2size=(12,2), fontsize=14,
    msize=10, showMLnoise=false, showresids=false, nresidbins=11)

    @assert burnin != 0
    M, = get_posterior(opt, fdataname=fdataname)
    from, = divrem(burnin, opt.save_freq)
    M = M[from+1:decimate:end]
    cols = 2
    showMLnoise && (cols=3; figsize=(15,7))
    fig, ax = plt[:subplots](1,cols, sharex=true, figsize=figsize)
    σ = zeros(size(d))
    cmap = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
    f3, ax3 = plt[:subplots](2,1, sharex=true, figsize=figsize)
    if showresids
        fig2, ax2 = plt[:subplots](1,length(opt_EM.freqs), figsize=fig2size, sharex=true, sharey=true)
        residedges = LinRange(-5.5,5.5,nresidbins)
        residbins = 0.5*(residedges[1:end-1]+residedges[2:end])
        histresid = zeros(length(residbins))
    end
    for (i, m) in enumerate(M)
        @assert length(opt_EM.rho) == length(m) + opt_EM.nfixed
        opt_EM.rho[opt_EM.nfixed+1:end] = 10.0.^m
        Er = CSEM1Dkernels.getCSEM1DanisoHED(opt_EM.freqs,
        opt_EM.rRx, opt_EM.zRx, opt_EM.zTx, opt_EM.z, [opt_EM.rho opt_EM.rho],
        0, RxAzim = opt_EM.RxAzim, TxDip = opt_EM.TxDip)

        for lfreq = 1:length(opt_EM.freqs)
            nr = sum(.!isnan.(d[:,lfreq]))
            select = .!isnan.(d[:,lfreq])
            r = Er[select,lfreq] - d[select,lfreq]
            σ[:,lfreq] = sqrt(0.5/nr)*norm(r./abs.(d[select,lfreq]), 2)*abs.(d[:,lfreq])
            ax3[1][:plot](opt_EM.rRx[select], real(r./σ[select,lfreq]), color=cmap[lfreq], alpha=alpha, ".", markersize=msize)
            ax3[2][:plot](opt_EM.rRx[select], imag(r./σ[select,lfreq]), color=cmap[lfreq], alpha=alpha, ".", markersize=msize)
            if showresids
                histresid[:] += fit(Histogram, real(r./σ[select,lfreq]), residedges).weights +
                                fit(Histogram, imag(r./σ[select,lfreq]), residedges).weights
                # histi[:] += fit(Histogram, imag(r./σ[select,lfreq]), residedges).weights
            end
            ax[1][:plot](opt_EM.rRx[select], abs.(Er[select,lfreq]), color=cmap[lfreq], alpha=alpha)
            ax[2][:plot](opt_EM.rRx[select], unwrap(angle.(Er[select,lfreq])), color=cmap[lfreq], alpha=alpha)
            if i == length(M)
                 ax[1][:errorbar](opt_EM.rRx[select], abs.(d[select,lfreq]), yerr = 2*sfactor*opt_EM.sd[select,lfreq],
                 linestyle="none", elinewidth=0, capsize=3, label=string(opt_EM.freqs[lfreq])*" Hz")
                ax[2][:errorbar](opt_EM.rRx[select], unwrap(angle.(d[select,lfreq])), yerr = 2*sfactor*opt_EM.sd[select,lfreq]./abs.(d[select,lfreq]),
                linestyle="none", elinewidth=0, capsize=3)
                showMLnoise && ax[3][:plot](opt_EM.rRx[select], sfactor*opt_EM.sd[select,lfreq], "x", markersize=msize)

                if showresids
                    histresid = histresid/sum(histresid)/(residbins[2]-residbins[1])
                    # histi = histi/sum(histi)/(residbins[2]-residbins[1])
                    ax2[lfreq][:plot](residbins, histresid, "-r")
                    # ax2[lfreq][:plot](residbins, histi, "-b")
                    ax2[lfreq][:plot](-6:.1:6, 1/sqrt(2*pi)*exp.(-0.5*(-6:.1:6).^2),"--k")
                    ax2[lfreq][:set_title]("$(opt_EM.freqs[lfreq]) Hz")
                    ax2[lfreq][:set_xlabel]("residuals")
                    ax2[lfreq][:set_ylim](0,0.55)
                    figure(fig2[:number])
                    MCMC_Driver.nicenup(gcf(),fsize=fontsize)
                end
            end
            showMLnoise && ax[3][:plot](opt_EM.rRx[select], σ[select,lfreq], ".", color=cmap[lfreq], markersize=msize, alpha=alpha)
        end
    end
    ax[1][:set_yscale]("log")
    ax[1][:grid]()
    ax[1][:set_xlabel]("offset m")
    ax[1][:set_ylabel]("Amplitude "* L"\frac{V}{Am^2}")
    ax[2][:grid]()
    ax[2][:set_xlabel]("offset m")
    ax[2][:set_ylabel]("Phase rad")
    if showMLnoise
        ax[3][:set_xlabel]("offset m")
        ax[3][:set_ylabel](L"\^{\sigma_{|amp|}}")
        ax[3][:grid]()
        ax[3][:set_yscale]("log")
        # @info("here23")
        # ax[3][:get_shared_y_axes]()[:join](ax[1], ax[3])
        # ax[3][:autoscale]()
    end

    figure(fig[:number])
    MCMC_Driver.nicenup(gcf(),fsize=fontsize)
    savefig(fdataname*"_post_fields.png", dpi=300)
    figure(fig2[:number])
    savefig(fdataname*"_post_resids.png", dpi=300)
    nothing
end

end
