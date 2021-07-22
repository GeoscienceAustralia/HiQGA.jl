module SurfaceRegression
import ..AbstractOperator.get_misfit
using ..AbstractOperator, ..CommonToAll
using PyPlot, LinearAlgebra, StatsBase

import ..Model, ..Options, ..OptionsStat, ..OptionsNonstat

export Surface, SurfaceWithDifferentData

mutable struct Surface <: Operator2D
    d      :: Array{Float64}
    useML  :: Bool
    σ      :: Array{Float64}
    r      :: Array{Float64}
    select :: Array{Bool} 
end

function Surface(d::Array{Float64} ;useML=false, σ=1.0)
    if length(d) != length(σ)
        @assert length(σ) == 1
        σ = σ[1]*ones(length(d))    
    end
    r = d[.!isnan.(d)]
    select = .!isnan.(d)
    Surface(d, useML, σ, r, select)
end

function get_misfit(m::Model, opt::Options, surface::Surface)
    chi2by2 = 0.0
    if !opt.debug
        d, σ, select, r  = surface.d, surface.σ, surface.select, surface.r
        r .= (m.fstar[vec(select)] .- d[select])./σ[select]
        chi2 = r'*r
        if surface.useML
            n = length(r)
            chi2by2 = 0.5*n*log(chi2)
        else
            chi2by2 = 0.5*chi2
        end
    end
    return chi2by2
end

mutable struct SurfaceWithDifferentData <: Operator2D
    surfaces :: Array{Surface, 1}
end    

function get_misfit(m::Model, opt::Options, surface::SurfaceWithDifferentData)
    chi2by2 = 0.0
    for s in surface.surfaces
        chi2by2 += get_misfit(m, opt, s)
    end    
    return chi2by2
end

function slice_image_posterior( M::AbstractArray, opt::Options, roworcol::Symbol;
                                rowcolnum = 1,
                                nbins = 50,
                                rhomin=Inf,
                                rhomax=-Inf,
                                qp1=0.05,
                                qp2=0.95,
                                pdfnormalize=false,
                                temperaturenum=1)
    @assert temperaturenum == 1
    Mslices = Array{Array{Float64}, 1}(undef, length(M))
    x = unique(opt.xall[1,:])
    y = unique(opt.xall[2,:])
    for i in 1:length(M)
        m = reshape(M[i],length(y), length(x))
        if roworcol == :row
            Mslices[i] = m[rowcolnum,:]
        else
            Mslices[i] = m[:,rowcolnum]
        end
    end
    himage, edges, CI, = gethimage(Mslices, opt, burninfrac=0.0, temperaturenum=temperaturenum,
                nbins=nbins, rhomin=rhomin, rhomax=rhomax, qp1=qp1, qp2=qp2,
                pdfnormalize=pdfnormalize)
end

function get_image_quantile(M::AbstractArray, q=0.5)
    n = length(M[1])
    inrows = size(M[1], 2) == 1
    if inrows 
        mall = hcat(M...)
    else 
        mall = vcat(M...)
    end        
    quantM = zeros(n)
    for i = 1:n
        quantM[i] = inrows ? quantile(mall[i,:], q) : quantile(mall[:,i], q)
    end
    quantM
end

function plot_surface_posterior(optns::OptionsNonstat, opts::OptionsStat;
                        rownum = 10,
                        colnum = 10,
                        temperaturenum = 1,
                        nbins = 50,
                        burninfrac=0.5,
                        qp1=0.05,
                        qp2=0.95,
                        cmappdf = "magma",
                        cmapmean = "bone",
                        figsizerows=[11,6],
                        figsizecols=[7.5,9.5],
                        pdfnormalize=false,
                        fsize=10,
                        domean = false,
                        getquant = 0.5,
                        concisefigsize=(10,5),
                        property="pixel value",
                        property_units="")
    @assert temperaturenum == 1
    @assert 0.0 < getquant < 1.0
    M = assembleTat1(optns, :fstar, burninfrac=burninfrac, temperaturenum=temperaturenum)
    if domean
        mns = mean(M)
    else
        mns = get_image_quantile(M, getquant)
    end
    x = vec(unique(optns.xall[1,:]))
    delx = diff(x)[1]
    xmesh = x[1]-delx/2:delx:x[end]+delx/2
    y = vec(unique(optns.xall[2,:]))
    dely = diff(y)[1]
    ymesh = y[1]-dely/2:dely:y[end]+dely/2
    mns = reshape(mns, length(y), length(x))
    CIimage = reshape(getCIimage(M, qp1, qp2), size(mns)) 
    himage_r_ns, edges_r_ns, CI_r_ns = slice_image_posterior(M, optns, :row, rowcolnum=rownum, nbins = nbins, qp1=qp1, qp2=qp2,
                                pdfnormalize=pdfnormalize, temperaturenum=temperaturenum)
    himage_c_ns, edges_c_ns, CI_c_ns = slice_image_posterior(M, optns, :col, rowcolnum=colnum, nbins = nbins, qp1=qp1, qp2=qp2,
                                pdfnormalize=pdfnormalize, temperaturenum=temperaturenum)
    M = assembleTat1(opts, :fstar, burninfrac=burninfrac, temperaturenum=temperaturenum)
    if domean
        ms = 0.5*log10.(mean(M))
    else
        a, b = [m[1,:] for m in M], [m[2,:] for m in M]
        ms = 0.5*log10.([get_image_quantile(a, getquant) get_image_quantile(b, getquant)])'
    end
    msx = reshape(ms[1,:], length(y), length(x))
    msy = reshape(ms[2,:], length(y), length(x))
    Mx, My = Array{Array{Float64}, 1}(undef, length(M)), Array{Array{Float64}, 1}(undef, length(M))
    for i = 1:length(M)
        Mx[i] = M[i][1,:]
        My[i] = M[i][2,:]
    end
    himage_x_c, edges_x_c, CI_x_c = slice_image_posterior(Mx, opts, :col, rowcolnum=colnum, nbins = nbins, qp1=qp1, qp2=qp2,
                                pdfnormalize=pdfnormalize, temperaturenum=temperaturenum)
    himage_y_c, edges_y_c, CI_y_c = slice_image_posterior(My, opts, :col, rowcolnum=colnum, nbins = nbins, qp1=qp1, qp2=qp2,
                                pdfnormalize=pdfnormalize, temperaturenum=temperaturenum)
    himage_x_r, edges_x_r, CI_x_r = slice_image_posterior(Mx, opts, :row, rowcolnum=rownum, nbins = nbins, qp1=qp1, qp2=qp2,
                                pdfnormalize=pdfnormalize, temperaturenum=temperaturenum)
    himage_y_r, edges_y_r, CI_y_r = slice_image_posterior(My, opts, :row, rowcolnum=rownum, nbins = nbins, qp1=qp1, qp2=qp2,
                                pdfnormalize=pdfnormalize, temperaturenum=temperaturenum)
    f = figure(figsize=figsizerows)
    s1 = subplot(231)
    im1 = s1.imshow(mns, extent=[xmesh[1], xmesh[end], ymesh[1], ymesh[end]], origin="lower", cmap=cmapmean)
    s1.plot([x[1], x[end]], [y[rownum], y[rownum]], "--k", alpha=0.5)
    central_tendency = domean ? ("Mean") : ("Percentile "*"$(round(Int, getquant*100))")
    s1.set_title(central_tendency*" "*property)
    s1.set_ylabel("distance y")
    s1.set_xlabel("distance x")
    cb1 = colorbar(im1, ax=s1)
    s2 = subplot(232, sharex=s1, sharey=s1)
    im2 = s2.imshow(msx, extent=[xmesh[1], xmesh[end], ymesh[1], ymesh[end]], origin="lower", cmap=cmapmean)
    s2.plot([x[1], x[end]], [y[rownum], y[rownum]], "--k", alpha=0.5)
    s2.set_title(central_tendency*" "*L"\log_{10}λ_x")
    s2.set_ylabel("distance y")
    s2.set_xlabel("distance x")
    cb2 = colorbar(im2, ax=s2)
    s3 = subplot(233, sharex=s1, sharey=s1)
    im3 = s3.imshow(msy, extent=[xmesh[1], xmesh[end], ymesh[1], ymesh[end]], origin="lower", cmap=cmapmean)
    s3.plot([x[1], x[end]], [y[rownum], y[rownum]], "--k", alpha=0.5)
    s3.set_title(central_tendency*" "*L"\log_{10}λ_y")
    s3.set_ylabel("distance y")
    s3.set_xlabel("distance x")
    cb3 = colorbar(im3, ax=s3)
    s4 = subplot(234, sharex=s1)
    im4 = s4.pcolormesh(xmesh[:], edges_r_ns[:], himage_r_ns', cmap=cmappdf)
    s4.plot(x, mns[rownum,:][:],"--w", alpha=0.5)
    s4.plot(x, CI_r_ns,"--c", alpha=0.5)
    s4.set_title(property*" PDF")
    s4.set_ylabel(property)
    s4.set_xlabel("distance x")
    cb4 = colorbar(im4, ax=s4)
    cb4.ax.set_title("pdf")
    s5 = subplot(235, sharex=s1)
    im5 = s5.pcolormesh(xmesh[:], edges_x_r[:], himage_x_r', cmap=cmappdf)
    cb5 = colorbar(im5, ax=s5)
    cb5.ax.set_title("pdf")
    s5.plot(x, msx[rownum,:][:],"--w", alpha=0.5)
    s5.plot(x, CI_x_r,"--c", alpha=0.5)
    s5.set_title(L"\log_{10}λ_x"*" PDF")
    s5.set_ylabel(L"\log_{10}λ_x")
    s5.set_xlabel("distance x")
    s6 = subplot(236, sharex=s1)
    im6 = s6.pcolormesh(xmesh[:], edges_y_r[:], himage_y_r', cmap=cmappdf)
    cb6 = colorbar(im6, ax=s6)
    cb6.ax.set_title("pdf")
    s6.plot(x, msy[rownum,:][:],"--w", alpha=0.5)
    s6.plot(x, CI_y_r,"--c", alpha=0.5)
    s6.set_title(L"\log_{10}λ_y"*" PDF")
    s6.set_ylabel(L"\log_{10}λ_y")
    s6.set_xlabel("distance x")
    nicenup(f, fsize=fsize)

    f = figure(figsize=figsizecols)
    s1 = subplot(321)
    im1 = s1.imshow(mns, extent=[xmesh[1], xmesh[end], ymesh[1], ymesh[end]], origin="lower", cmap=cmapmean)
    s1.plot([y[colnum], y[colnum]], [y[1], y[end]], "--k", alpha=0.5)
    s1.set_title(central_tendency*" "*property)
    s1.set_ylabel("distance y")
    s1.set_xlabel("distance x")
    # s1.set_aspect(aspect)
    cb1 = colorbar(im1, ax=s1)
    s2 = subplot(323, sharex=s1, sharey=s1)
    im2 = s2.imshow(msx, extent=[xmesh[1], xmesh[end], ymesh[1], ymesh[end]], origin="lower", cmap=cmapmean)
    s2.plot([y[colnum], y[colnum]], [y[1], y[end]], "--k", alpha=0.5)
    s2.set_title(central_tendency*" "*L"\log_{10}λ_x")
    s2.set_ylabel("distance y")
    s2.set_xlabel("distance x")
    # s2.set_aspect(aspect)
    s3 = subplot(325, sharex=s1, sharey=s1)
    cb2 = colorbar(im2, ax=s2)
    im3 = s3.imshow(msy, extent=[xmesh[1], xmesh[end], ymesh[1], ymesh[end]], origin="lower", cmap=cmapmean)
    s3.plot([y[colnum], y[colnum]], [y[1], y[end]], "--k", alpha=0.5)
    # s3.set_aspect(aspect)
    s3.set_title(central_tendency*" "*L"\log_{10}λ_y")
    s3.set_ylabel("distance y")
    s3.set_xlabel("distance x")
    cb3 = colorbar(im3, ax=s3)
    s4 = subplot(322, sharey=s1)
    im4 = s4.pcolormesh(edges_c_ns[:], ymesh[:], himage_c_ns, cmap=cmappdf)
    s4.plot(mns[:,colnum][:], y, "--w", alpha=0.5)
    s4.plot(CI_c_ns, y, "--c", alpha=0.5)
    s4.set_title(property*" PDF")
    s4.set_xlabel(property)
    s4.set_ylabel("distance y")
    cb4 = colorbar(im4, ax=s4)
    cb4.ax.set_title("pdf")
    s5 = subplot(324, sharey=s1)
    im5 = s5.pcolormesh(edges_x_c[:], ymesh[:], himage_x_c, cmap=cmappdf)
    cb5 = colorbar(im5, ax=s5)
    cb5.ax.set_title("pdf")
    s5.plot(msx[:,colnum][:], y, "--w", alpha=0.5)
    s5.plot(CI_x_c, y, "--c", alpha=0.5)
    s5.set_title(L"\log_{10}λ_x"*" PDF")
    s5.set_xlabel(L"\log_{10}λ_x")
    s5.set_ylabel("distance y")
    s6 = subplot(326, sharey=s1)
    im6 = s6.pcolormesh(edges_y_c[:], ymesh[:], himage_y_c, cmap=cmappdf)
    cb6 = colorbar(im6, ax=s6)
    cb6.ax.set_title("pdf")
    s6.plot(msy[:,colnum][:], y, "--w", alpha=0.5)
    s6.plot(CI_y_c, y, "--c", alpha=0.5)
    s6.set_title(L"\log_{10}λ_y"*" PDF")
    s6.set_xlabel(L"\log_{10}λ_y")
    s6.set_ylabel("distance y")
    nicenup(f, fsize=fsize)

    f, ax = plt.subplots(1, 2, sharex=true, sharey=true, figsize=concisefigsize)
    im1 = ax[1].imshow(mns, extent=[xmesh[1], xmesh[end], ymesh[1], ymesh[end]], origin="lower", cmap=cmapmean)
    central_tendency = domean ? ("Mean") : ("Percentile "*"$(round(Int, getquant*100))")
    ax[1].set_title(central_tendency*" "*property)
    ax[1].set_ylabel("distance y")
    ax[1].set_xlabel("distance x")
    cb1 = colorbar(im1, ax=ax[1])
    cb1.ax.set_xlabel(property_units)
    im2 = ax[2].imshow(CIimage, extent=[xmesh[1], xmesh[end], ymesh[1], ymesh[end]], origin="lower", cmap=cmappdf)
    ax[2].set_title("Credible interval: $(round(Int, 100*qp2)) - $(round(Int, 100*qp1))")
    ax[2].set_xlabel("distance x")
    cb2 = colorbar(im2, ax=ax[2])
    cb2.ax.set_xlabel(property_units)
end

function getCIimage(marray, qp1, qp2)
    q1 = get_image_quantile(marray, qp1)
    q2 = get_image_quantile(marray, qp2)
    return q2-q1
end   


function plot_surface_posterior(opt::Options;
                        rownum = 10,
                        colnum = 10,
                        temperaturenum = 1,
                        nbins = 50,
                        burninfrac=0.5,
                        qp1=0.05,
                        qp2=0.95,
                        cmappdf = "magma",
                        cmapmean = "bone",
                        figsize=[7.5,5.9],
                        pdfnormalize=false,
                        fsize=10,
                        aspect=1.0,
                        domean = false,
                        getquant = 0.5,
                        property = "pixel value",
                        property_units="",
                        concisefigsize=(10,5))
    @assert temperaturenum == 1
    @assert 0.0 < getquant < 1.0
    if typeof(opt) == OptionsStat
        @assert opt.updatenonstat == false
        @assert opt.needλ²fromlog == false
    end
    x = vec(unique(opt.xall[1,:]))
    delx = diff(x)[1]
    xmesh = x[1]-delx/2:delx:x[end]+delx/2
    y = vec(unique(opt.xall[2,:]))
    dely = diff(y)[1]
    ymesh = y[1]-dely/2:dely:y[end]+dely/2
    M = assembleTat1(opt, :fstar, burninfrac=burninfrac, temperaturenum=temperaturenum)
    if domean
        m = mean(M)
    else
        a = [m' for m in M]
        m = get_image_quantile(a, getquant)
    end
    m = reshape(m, length(y), length(x))
    CIimage = reshape(getCIimage(M, qp1, qp2), size(m)) 
    himage_r, edges_r, CI_r = slice_image_posterior(M, opt, :row, rowcolnum=rownum, nbins = nbins, qp1=qp1, qp2=qp2,
                                pdfnormalize=pdfnormalize, temperaturenum=temperaturenum)
    himage_c, edges_c, CI_c = slice_image_posterior(M, opt, :col, rowcolnum=colnum, nbins = nbins, qp1=qp1, qp2=qp2,
                                pdfnormalize=pdfnormalize, temperaturenum=temperaturenum)
    f = figure(figsize=figsize)
    s1 = subplot(223)
    central_tendency = domean ? ("Mean") : ("Percentile "*"$(round(Int, getquant*100))")
    im1 = s1.imshow(m, extent=[xmesh[1], xmesh[end], ymesh[1], ymesh[end]], origin="lower", cmap=cmapmean)
    s1.plot([x[1], x[end]], [y[rownum], y[rownum]], "--k", alpha=0.5)
    s1.plot([y[colnum], y[colnum]], [y[1], y[end]], "--k", alpha=0.5)
    s1.set_title(central_tendency*" "*property)
    s1.set_ylabel("distance y")
    # s1.set_aspect(aspect)
    s1.set_xlabel("distance x")
    cb1 = colorbar(im1, ax=s1)
    s2 = subplot(221, sharex=s1)
    im2 = s2.pcolormesh(xmesh[:], edges_r[:], himage_r', cmap=cmappdf)
    s2.plot(x, m[rownum,:][:],"--w", alpha=0.5)
    s2.plot(x, CI_r, "--c", alpha=0.5)
    s2.set_title(property*" PDF")
    s2.set_ylabel(property)
    s2.set_xlabel("distance x")
    cb2 = colorbar(im2, ax=s2)
    cb2.ax.set_title("pdf")
    s3 = subplot(224, sharey=s1)
    im3 = s3.pcolormesh(edges_c[:], ymesh[:], himage_c, cmap=cmappdf)
    s3.plot(m[:,colnum][:], y, "--w", alpha=0.5)
    s3.plot(CI_c, y, "--c", alpha=0.5)
    s3.set_title(property*" PDF")
    s3.set_xlabel(property)
    s3.set_ylabel("distance y")
    cb3 = colorbar(im3, ax=s3)
    cb3.ax.set_title("pdf")
    nicenup(f, fsize=fsize)

    f, ax = plt.subplots(1, 2, sharex=true, sharey=true, figsize=concisefigsize)
    im1 = ax[1].imshow(m, extent=[xmesh[1], xmesh[end], ymesh[1], ymesh[end]], origin="lower", cmap=cmapmean)
    central_tendency = domean ? ("Mean") : ("Percentile "*"$(round(Int, getquant*100))")
    ax[1].set_title(central_tendency*" "*property)
    ax[1].set_ylabel("distance y")
    ax[1].set_xlabel("distance x")
    cb1 = colorbar(im1, ax=ax[1])
    cb1.ax.set_xlabel(property_units)
    im2 = ax[2].imshow(CIimage, extent=[xmesh[1], xmesh[end], ymesh[1], ymesh[end]], origin="lower", cmap=cmappdf)
    ax[2].set_title("Credible interval: $(round(Int, 100*qp2)) - $(round(Int, 100*qp1))")
    ax[2].set_xlabel("distance x")
    cb2 = colorbar(im2, ax=ax[2])
    cb2.ax.set_xlabel(property_units)
end

end
