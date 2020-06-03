module LineRegression
import AbstractOperator.get_misfit
import AbstractOperator.Operator
using TransD_GP, PyPlot, StatsBase, Statistics, LinearAlgebra, CommonToAll

export Line, plot_posterior, make1Dhist, make1Dhists

mutable struct Line<:Operator
    d     :: Array{Float64}
    useML :: Bool
    σ     :: Float64
end

function Line(d::Array{Float64, 1} ;useML=false, σ=1.0)
    Line(d, useML, σ)
end

function plot_posterior(L::Line,
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

function plot_posterior(L::Line,
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
    im1 = ax.imshow(himage, extent=[edges[1],edges[end],xall[end],xall[1]], aspect="auto", cmap=cmappdf)
    ax.grid()
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

function get_misfit(m::TransD_GP.Model, opt::TransD_GP.Options, line::Line)
    chi2by2 = 0.0
    if !opt.debug
        d = line.d
        select = .!isnan.(d[:])
        r = m.fstar[:][select] - d[select]
        if line.useML
            N = sum(select)
            chi2by2 = 0.5*N*log(norm(r)^2)
        else
            chi2by2 = r'*r/(2line.σ^2)
        end
    end
    return chi2by2
end

end
