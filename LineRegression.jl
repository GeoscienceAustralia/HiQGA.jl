using TransD_GP, PyPlot, StatsBase, Statistics, LinearAlgebra

mutable struct Line<:Operator
    n::Nothing
end

# function get_posterior(opt_in::TransD_GP.Options;
#     fdataname="",
#     burnin=300000,
#     nbins=50,
#     rho=[0.],
#     alpha=0.025,
#     alphacolor = "blue",
#     alphalwidth = 0.3,
#     qp1=0.01,
#     qp2=0.99,
#     figsize = (10,6),
#     fontsize = 16,
#     convfigsize = (8,6),
#     kfigsize = (9,6),
#     convfontsize = 14,
#     cmapcdf = "RdYlBu_r",
#     cmappdf = "inferno_r",
#     showtitle=false,
#     titlestr="",
#     topadjust=0.9,
#     fixviewextent=[0. 0. 0. 0.],
#     rhobgcolor="yellow",
#     rhofgcolor="k",
#     nlldelta = 100,
#     nxbins = 20,
#     nftbins = 20,
#     nxftfontsize = 16,
#     pdfnormalize=false)
#     @assert length(fixviewextent) == 4
#
#     M, n, x_ft, iter, misfit, T0loc = get_posterior(opt_in, fdataname=fdataname)
#
#     from = findfirst(iter.>burnin)
#     f1 = figure(figsize=convfigsize)
#     s1 = subplot(311)
#     plot(iter,n)
#     s1.set_ylim(opt_in.nmin, opt_in.nmax)
#     ylabel("# train")
#     s2 = subplot(312, sharex=s1)
#     plot(iter, misfit)
#     ylabel("-log L")
#     s3 = subplot(313, sharex=s1)
#     plot(iter, T0loc)
#     s3.set_yticks(collect(minimum(T0loc):2:maximum(T0loc)))
#     s2.set_ylim(mean(misfit[from:end])-nlldelta, mean(misfit[from:end])+nlldelta)
#     s3.set_xlim(burnin+1, iter[end])
#     ylabel("T=1 location")
#     xlabel("iterations")
#     MCMC_Driver.nicenup(gcf(), fsize=convfontsize)
#     if showtitle
#         titlestr == "" && (titlestr = fdataname)
#         f1.suptitle(titlestr, fontsize=convfontsize)
#         f1.subplots_adjust(top=topadjust)
#     end
#
#     xall = opt_in.xall
#     f,ax = plt.subplots(1,2, sharex=true, sharey=true, figsize=figsize)
#     rhomin, rhomax = 0.0, 0.0
#     for (i,mm) in enumerate(M)
#         if i >= from
#             # ax[1].step(mm, xall[:], linewidth=alphalwidth, color=alphacolor, alpha=alpha)
#             #plot(x_ft[i][1:n[i],2], x_ft[i][1:n[i],1], ".", markersize=20)
#             # title("model $(iter[i])")
#         end
#         rhomin_mm = min(mm...)
#         rhomax_mm = max(mm...)
#         rhomin_mm < rhomin && (rhomin = rhomin_mm)
#         rhomax_mm > rhomax && (rhomax = rhomax_mm)
#     end
#
#     # if rho != [0.]
#     #     ax[1].step(log10.(rho[2:end]), opt_EM.z[2:end], linewidth=3, color=rhobgcolor)
#     #     ax[1].step(log10.(rho[2:end]), opt_EM.z[2:end], linewidth=2, color=rhofgcolor)
#     # end
#     # ax[1][:grid]()
#     # ax[1].invert_yaxis()
#
#     edges = LinRange(rhomin, rhomax, nbins+1)
#     himage = zeros(Float64, length(M[1]), nbins)
#     cumhimage = zeros(Float64, length(M[1]), nbins)
#     CI = zeros(Float64, length(M[1]), 2)
#     z = opt_EM.z
#     for ilayer=1:length(M[1])
#         himage[ilayer,:] = (fit(Histogram, [m[ilayer] for m in M[from:end]], edges).weights)
#         himage[ilayer,:] = himage[ilayer,:]/sum(himage[ilayer,:])/(diff(edges)[1])
#         cumhimage[ilayer,:] = cumsum(himage[ilayer,:])/sum(himage[ilayer,:])
#         pdfnormalize && (himage[ilayer,:] = himage[ilayer,:]/maximum(himage[ilayer,:]))
#         CI[ilayer,:] = [quantile([m[ilayer] for m in M[from:end]],(qp1, qp2))...]
#     end
#     im1 = ax[1][:imshow](himage, extent=[edges[1],edges[end],xall[end],xall[1]], aspect="auto",
#             cmap=cmappdf)
#     cb1 = colorbar(im1, ax=ax[1])
#     cb1[:ax].set_xlabel("pdf")
#     ax[1][:grid]()
#     im2 = ax[2][:imshow](cumhimage, extent=[edges[1],edges[end],xall[end],xall[1]], aspect="auto",
#             cmap=cmapcdf)
#     ax[2][:grid]()
#     cb2 = colorbar(im2, ax=ax[2])
#     cb2[:ax].set_xlabel("cdf")
#     if rho != [0.]
#         ax[1].step(log10.(rho[2:end]), opt_EM.z[2:end], linewidth=3, color=rhobgcolor)
#         ax[1].step(log10.(rho[2:end]), opt_EM.z[2:end], linewidth=2, color=rhofgcolor)
#         ax[2].step(log10.(rho[2:end]), opt_EM.z[2:end], linewidth=3, color=rhobgcolor)
#         ax[2].step(log10.(rho[2:end]), opt_EM.z[2:end], linewidth=2, color=rhofgcolor)
#     end
#     ax[1].plot(CI, xall[:], linewidth=2, color="g")
#     ax[1].plot(CI, xall[:], linewidth=2, color="k", linestyle="--")
#     ax[2].plot(CI, xall[:], linewidth=2, color="g")
#     ax[2].plot(CI, xall[:], linewidth=2, color="k", linestyle="--")
#     # ax[2].plot(CI, xall[:], "--k")
#     ylim(opt_EM.zRx, xall[end])
#     gca().invert_yaxis()
#     ax[1].set_xlabel(L"\log_{10} \rho")
#     ax[1].set_ylabel("depth (m)")
#     ax[2].set_xlabel(L"\log_{10} \rho")
#     if !(fixviewextent == [0. 0. 0. 0.])
#         ax[1].set_xlim(fixviewextent[1:2])
#         ax[1].set_ylim(fixviewextent[3:4])
#         ax[1].invert_yaxis()
#     end
#     MCMC_Driver.nicenup(f, fsize=fontsize)
#     if showtitle
#         titlestr == "" && (titlestr = fdataname)
#         f.suptitle(titlestr, fontsize=fontsize)
#         f.subplots_adjust(top=topadjust)
#     end
#
#     x, ft = zeros(sum(n[from:end])), zeros(sum(n[from:end]))
#     nlast = 0
#     for (i, xt) in enumerate(x_ft[from:end])
#         x[nlast+1:nlast+n[from-1+i]]  = xt[1:n[from-1+i],1]
#         ft[nlast+1:nlast+n[from-1+i]] = xt[1:n[from-1+i],2]
#         nlast = nlast+n[from-1+i]
#     end
#     f2, ax2 = plt.subplots(2,2, figsize=kfigsize)
#     edgesx = LinRange(opt_in.xbounds[1], opt_in.xbounds[2], nxbins+1)
#     edgesrho = LinRange(opt_in.fbounds[1], opt_in.fbounds[2], nftbins+1)
#     h = fit(Histogram, (x, ft), (edgesx, edgesrho)).weights
#     im3 = ax2[1][:imshow](h, extent=[edgesrho[1],edgesrho[end],edgesx[end],edgesx[1]], aspect="auto", vmin=0.0, vmax=2*max(h...))
#     ax2[1].set_ylabel("depth m")
#     ax2[1].set_xlabel(L"\log_{10}\rho")
#     rhist = sum(h,dims=1)[:]
#     rhist = rhist/sum(rhist)/diff(edgesrho)[1]
#     ax2[2][:bar](0.5*(edgesrho[2:end]+edgesrho[1:end-1]), rhist, width=diff(edgesrho)[1])
#     ax2[2].plot([edgesrho[1], edgesrho[end]], 1/(opt_in.fbounds[end]-opt_in.fbounds[1])*[1,1],"--k")
#     ax2[2].set_ylim(0, 3*maximum(rhist))
#     ax2[2].set_xlabel(L"\log_{10}\rho")
#     ax2[2].set_ylabel("pdf")
#     xhist = sum(h,dims=2)[:]
#     xhist = xhist/sum(xhist)/diff(edgesx)[1]
#     ax2[3][:barh](0.5*(edgesx[2:end]+edgesx[1:end-1]), xhist, height=diff(edgesx)[1])
#     ax2[3].plot(1/(opt_in.xbounds[end]-opt_in.xbounds[1])*[1,1], [edgesx[1], edgesx[end]],"--k")
#     ax2[3].set_xlim(0, 3*maximum(xhist))
#     ax2[3].set_ylim(ax2[1][:get_ylim]())
#     ax2[3].set_ylabel("depth m")
#     ax2[3].set_xlabel("pdf")
#     ax2[2][:get_shared_y_axes]()[:join](ax2[1], ax2[3])
#     ax2[2][:get_shared_x_axes]()[:join](ax2[1], ax2[2])
#     ax2[2].set_xlim(ax2[1][:get_xlim]())
#     k = fit(Histogram, n[from:end], (opt_in.nmin-0.5:opt_in.nmax+0.5)).weights
#     k = k/sum(k)
#     @show mean(k)
#     ax2[4][:bar](opt_in.nmin:opt_in.nmax, k, width=1)
#     ax2[4].plot([opt_in.nmin, opt_in.nmax],1/(opt_in.nmax-opt_in.nmin+1)*[1,1],"--k")
#     ax2[4].set_xlabel("# training points")
#     ax2[4].set_ylim(0, 3*max(k...))
#     ax2[4].set_ylabel("pdf")
#     ax2[4].set_xlim(opt_in.nmin, opt_in.nmax)
#     MCMC_Driver.nicenup(gcf(), fsize=nxftfontsize)
#     savefig(fdataname*"_k.png", dpi=300)
#
#     figure(f1.number)
#     savefig(fdataname*"_conv.png", dpi=300)
#     figure(f.number)
#     savefig(fdataname*"_post.png", dpi=300)
#     return ax[1], im1, cb1, ax[2], im2, cb2
# end

function get_posterior(opt_in::TransD_GP.Options; fdataname="", decimate=1)

    M = TransD_GP.history(opt_in, stat=:fstar)
    n = TransD_GP.history(opt_in, stat=:nodes)
    x_ft = TransD_GP.history(opt_in, stat=:x_ftrain)
    iter = TransD_GP.history(opt_in, stat=:iter)

    return M[1:decimate:end], n[1:decimate:end],
    x_ft[1:decimate:end], iter[1:decimate:end]
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
