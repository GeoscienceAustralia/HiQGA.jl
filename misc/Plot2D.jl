module Plot2D
using TransD_GP, PyPlot, MCMC_Driver, JLD, StatsBase, Statistics
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

function plot_posterior(X::Array{Float64,1}, Y::Array{Float64,1}, Xall::Array{Float64, 2},
                        opt_in::TransD_GP.Options;
        burnin=50000,
        nbins=50,
        figsize=(8,8),
        fdataname="",
        fontsize=14)

    M, n, x_ft, iter, misfit, T0loc = get_posterior(opt_in, fdataname=fdataname)

    f2, ax2 = plt[:subplots](2,1, sharex=true, figsize=(6,4))
    ax2[1][:plot](iter,n)
    ax2[1][:grid]()
    ax2[2][:plot](iter, misfit)
    ax2[2][:grid]()
    ax2[2][:set_xlabel]("iterations")
    ax2[2][:set_ylabel]("-log L")
    ax2[1][:set_ylabel]("# training")
    MCMC_Driver.nicenup(gcf(), fsize=13)
    savefig(fdataname*"_conv.png", dpi=300)

    from = findfirst(iter.>burnin)

    xall = opt_in.xall
    sq = zeros(size(M[1]))
    mn = zeros(size(M[1]))
    for mm in M[from:end]
        sq += mm.^2
        mn += mm
    end

    N = length(M[from:end])
    mn = mn/N

    f3, ax3 = plt[:subplots](1,2,figsize=(10,5), sharex=true, sharey=true)
    nmodel = from
    im1 = ax3[1][:imshow](reshape(M[nmodel],length(Y), length(X)), extent=[X[1],X[end],Y[end],Y[1]])
    ax3[1][:scatter](x_ft[nmodel][1:n[nmodel],1], x_ft[nmodel][1:n[nmodel],2], s=20, color="black", alpha=0.5)
    cb1 = colorbar(im1, ax=ax3[1])
    im2 = ax3[2][:imshow](reshape(mn,length(Y), length(X)), extent=[X[1],X[end],Y[end],Y[1]])
    cb2 = colorbar(im2, ax=ax3[2])
    MCMC_Driver.nicenup(gcf(), fsize=14)
    savefig(fdataname*"_final_mean.png", dpi=300)

    vr = sq/N - (mn).^2

    x, y = zeros(sum(n[from:end])), zeros(sum(n[from:end]))
    nlast = 0
    for (i, xt) in enumerate(x_ft[from:end])
        x[nlast+1:nlast+n[from-1+i]]  = xt[1:n[from-1+i],1]
        y[nlast+1:nlast+n[from-1+i]] = xt[1:n[from-1+i],2]
        nlast = nlast+n[from-1+i]
    end

    f2, ax2 = plt[:subplots](2,2, figsize=figsize)
    edgesx = LinRange(opt_in.xbounds[1,1], opt_in.xbounds[1,2], nbins+1)
    edgesy = LinRange(opt_in.xbounds[2,1], opt_in.xbounds[2,2], nbins+1)
    h = fit(Histogram, (y,x), (edgesy, edgesx)).weights
    im3 = ax2[1][:imshow](h, extent=[edgesy[1],edgesy[end],edgesx[end],edgesx[1]], aspect="auto")
    # ax2[1][:set_ylabel]("m")
    # ax2[1][:set_xlabel]("m")
    yhist = sum(h,dims=1)[:]
    yhist = yhist/sum(yhist)/diff(edgesy)[1]
    ax2[2][:bar](0.5*(edgesy[2:end]+edgesy[1:end-1]), yhist, width=diff(edgesy)[1])
    ax2[2][:set_xlabel]("y")
    ax2[2][:set_ylabel]("pdf")
    xhist = sum(h,dims=2)[:]
    xhist = xhist/sum(xhist)/diff(edgesx)[1]
    ax2[3][:barh](0.5*(edgesx[2:end]+edgesx[1:end-1]), xhist, height=diff(edgesx)[1])
    ax2[3][:set_xlim](0, 3*maximum(xhist))
    ax2[3][:set_ylim](ax2[1][:get_ylim]())
    ax2[3][:set_ylabel]("x")
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
    ax2[4][:set_ylabel]("pdf")
    ax2[4][:set_xlim](opt_in.nmin, opt_in.nmax)
    MCMC_Driver.nicenup(gcf(), fsize=fontsize)
    savefig(fdataname*"_post_k.png", dpi=300)

    f,ax = plt[:subplots](1,2,figsize=(10,5), sharex=true, sharey=true)
    im1 = ax[1][:imshow](reshape(log10.(sqrt.(vr)),length(Y), length(X)), extent=[X[1],X[end],Y[end],Y[1]], cmap="gray_r")
    ax[1][:scatter](Xall[1,:], Xall[2,:], s=5, color="red", alpha=0.25)
    cb1 = colorbar(im1, ax=ax[1])
    im2 = ax[2][:imshow](log10.(h), extent=[edgesy[1],edgesy[end],edgesx[end],edgesx[1]], cmap="gray_r")
    cb2 = colorbar(im2, ax=ax[2])
    @info "lala6"
    # ax[1][:set_ylabel]("m")
    # ax[1][:set_xlabel]("m")
    # ax[2][:set_xlabel]("m")
    MCMC_Driver.nicenup(gcf(), fsize=14)
    savefig(fdataname*"_final_sd", dpi=300)

    return
end

end
