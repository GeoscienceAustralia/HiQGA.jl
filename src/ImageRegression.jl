module ImageRegression
import AbstractOperator.get_misfit
using AbstractOperator
using TransD_GP, PyPlot, LinearAlgebra, Statistics
using Random, Images, CommonToAll

export Img, get_image_data, calc_image_RMS, get_image_prediction_points,
        plot_image_data, plot_image_posterior, slice_image_posterior

mutable struct Img<:Operator2D

    filename          :: String
    x                 :: StepRangeLen
    y                 :: StepRangeLen
    fractrain         :: Float64
    dec               :: Int
    gausskernelwidth  :: Int
    f                 :: Array{Float64, 2}
    d                 :: Array{Float64, 2}
    useML             :: Bool
    σ                 :: Float64
end

function Img(;
             filename         = "",
             dx               = 0.01,
             fractrain        = 0.02,
             dec::Int         = 2,
             gausskernelwidth = 7,
             useML            = false,
             σ                = 1.0)

    @assert fractrain > 0 && fractrain < 1
    @assert !(filename == "")
    @assert dec > 1
    @assert σ > 0.0
    f = get_image(filename, gausskernelwidth, dec)
    x = 0:dx:dx*size(f,2)-1
    y = 0:dx:dx*size(f,1)-1
    Img(filename, x, y, fractrain, dec, gausskernelwidth, f, f, useML, σ)
end

function get_image(filename::String, gausskernelwidth::Int, dec)
    f = Gray.(load(filename))
    f = convert(Array{Float64, 2}, f)[1:dec:end,:1:dec:end]
    # convert image into something that looks like log resistivity
    f = -1 .+ 3*f
    # smooth image
    imfilter(f,Kernel.gaussian(gausskernelwidth))
end

function get_image_data(img::Img;
                           sdmaxfrac = 0.05,
                           rseed     = 12,
                           ybreak    = 12321.0,
                           takeevery = 4
                          )
    @assert ybreak != 12321.0
    @assert sdmaxfrac > 0 && sdmaxfrac < 1
    Random.seed!(rseed)
    f, x, y = img.f, img.x, img.y
    noisyd = NaN .+ zeros(Float64, size(f))
    ntrain = round(Int, img.fractrain*length(f))
    δtry = sdmaxfrac*max(f...)
    Xtrain = zeros(2,0)
    linidx = randperm(length(f))[1:ntrain]
    lgood = zeros(Int, 0)
    for (i,l) in enumerate(linidx)
        row, col = Tuple(CartesianIndices(f)[l])
        if img.y[row] > ybreak
            noisyd[row, col] = f[row, col] + δtry*randn()
            push!(lgood, l)
            Xtrain = hcat(Xtrain, [x[col]; y[row]])
        else
            if rem(row, takeevery) == 0 && rem(col, takeevery) == 0
                noisyd[row, col] = f[row, col] + δtry*randn()
                push!(lgood, l)
                Xtrain = hcat(Xtrain, [x[col]; y[row]])
            end
        end
    end
    img.d = noisyd
    δtry, noisyd[lgood], Xtrain
end

function get_image_prediction_points(img::Img)
    f, x, y = img.f, img.x, img.y
    Xall = zeros(2, length(f))
    for i in 1:length(f)
        yid, xid = Tuple(CartesianIndices(f)[i])
        Xall[:,i] = [x[xid]; y[yid]]
    end
    Xall
end

function plot_image_data(ftrain::Array{Float64, 1}, Xtrain::Array{Float64, 2},
                   img::Img;  s=3, fsize=12)
    f, x, y = img.f, img.x, img.y
    f1, ax1 = plt.subplots(1,2,figsize=(6.7,2.7), sharex=true, sharey=true)
    im1 = ax1[1].imshow(f, extent=[x[1],x[end],y[end],y[1]])
    cb1 = colorbar(im1, ax=ax1[1])
    ax1[2].imshow(f, extent=[x[1],x[end],y[end],y[1]], alpha=0.0)
    im2 = ax1[2].scatter(Xtrain[1,:], Xtrain[2,:], c=ftrain, s=s)
    cb2 = colorbar(im2)
    ax1[2].axis([x[1],x[end],y[end],y[1]])
    nicenup(gcf(), fsize=fsize)

end

function calc_image_RMS(img::Img)
    select = .!isnan.(img.d)
    r = (img.d[select] - img.f[select])/img.σ
    n = sum(select)
    @info "χ^2 error is $(r'*r) for $n points RMS: $(sqrt(r'*r/n))"
    nothing
end

function get_misfit(m::TransD_GP.Model, opt::TransD_GP.Options, img::Img)
    chi2by2 = 0.0
    if !opt.debug
        d = img.d
        select = .!isnan.(d[:])
        r = m.fstar[select] - d[select]
        if img.useML
            N = sum(select)
            chi2by2 = 0.5*N*log(norm(r)^2)
        else
            chi2by2 = r'*r/(2img.σ^2)
        end
    end
    return chi2by2
end

function slice_image_posterior( M::AbstractArray,
                                opt::TransD_GP.Options,
                                img::Img, roworcol::Symbol;
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
    for i in 1:length(M)
        m = reshape(M[i],length(img.y), length(img.x))
        if roworcol == :row
            Mslices[i] = m[rowcolnum,:]
        else
            Mslices[i] = m[:,rowcolnum]
        end
    end
    himage, edges, CI = gethimage(Mslices, opt, burninfrac=0.0, temperaturenum=temperaturenum,
                nbins=nbins, rhomin=rhomin, rhomax=rhomax, qp1=qp1, qp2=qp2,
                pdfnormalize=pdfnormalize)
end

function get_image_quantile(M::AbstractArray, q=0.5)
    mall = hcat(M...)
    quantM = zeros(length(M[1]))
    for i = 1:length(M[1])
        quantM[i] = quantile(mall[i,:], q)
    end
    quantM
end

function plot_image_posterior(optns::TransD_GP.OptionsNonstat, opts::TransD_GP.OptionsStat,
                        img::Img;
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
                        aspect=1.0,
                        domean = false,
                        getquant = 0.5)
    @assert temperaturenum == 1
    @assert 0.0 < getquant < 1.0
    M = assembleTat1(optns, :fstar, burninfrac=burninfrac, temperaturenum=temperaturenum)
    if domean
        mns = mean(M)
    else
        mns = get_image_quantile(M, getquant)
    end
    mns = reshape(mns, length(img.y), length(img.x))
    himage_r_ns, edges_r_ns, CI_r_ns = slice_image_posterior(M, optns, img, :row, rowcolnum=rownum, nbins = nbins, qp1=qp1, qp2=qp2,
                                pdfnormalize=pdfnormalize, temperaturenum=temperaturenum)
    himage_c_ns, edges_c_ns, CI_c_ns = slice_image_posterior(M, optns, img, :col, rowcolnum=colnum, nbins = nbins, qp1=qp1, qp2=qp2,
                                pdfnormalize=pdfnormalize, temperaturenum=temperaturenum)
    M = assembleTat1(opts, :fstar, burninfrac=burninfrac, temperaturenum=temperaturenum)
    if domean
        ms = 0.5*log10.(mean(M))
    else
        a, b = [m[1,:] for m in M], [m[2,:] for m in M]
        ms = 0.5*log10.([get_image_quantile(a, getquant) get_image_quantile(b, getquant)])'
    end
    msx = reshape(ms[1,:], length(img.y), length(img.x))
    msy = reshape(ms[2,:], length(img.y), length(img.x))
    Mx, My = Array{Array{Float64}, 1}(undef, length(M)), Array{Array{Float64}, 1}(undef, length(M))
    for i = 1:length(M)
        Mx[i] = M[i][1,:]
        My[i] = M[i][2,:]
    end
    himage_x_c, edges_x_c, CI_x_c = slice_image_posterior(Mx, opts, img, :col, rowcolnum=colnum, nbins = nbins, qp1=qp1, qp2=qp2,
                                pdfnormalize=pdfnormalize, temperaturenum=temperaturenum)
    himage_y_c, edges_y_c, CI_y_c = slice_image_posterior(My, opts, img, :col, rowcolnum=colnum, nbins = nbins, qp1=qp1, qp2=qp2,
                                pdfnormalize=pdfnormalize, temperaturenum=temperaturenum)
    himage_x_r, edges_x_r, CI_x_r = slice_image_posterior(Mx, opts, img, :row, rowcolnum=rownum, nbins = nbins, qp1=qp1, qp2=qp2,
                                pdfnormalize=pdfnormalize, temperaturenum=temperaturenum)
    himage_y_r, edges_y_r, CI_y_r = slice_image_posterior(My, opts, img, :row, rowcolnum=rownum, nbins = nbins, qp1=qp1, qp2=qp2,
                                pdfnormalize=pdfnormalize, temperaturenum=temperaturenum)
    f = figure(figsize=figsizerows)
    s1 = subplot(231)
    im1 = s1.pcolormesh(img.x, img.y, mns, cmap=cmapmean)
    s1.plot([img.x[1], img.x[end]], [img.y[rownum], img.y[rownum]], "--k", alpha=0.5)
    central_tendency = domean ? ("Mean") : ("Percentile "*"$(round(Int, getquant*100))")
    s1.set_title(central_tendency*" pixel value")
    s1.set_ylabel("distance y")
    # s1.set_aspect(aspect)
    s1.set_xlabel("distance x")
    cb1 = colorbar(im1, ax=s1)
    s2 = subplot(232, sharex=s1, sharey=s1)
    im2 = s2.pcolormesh(img.x, img.y, msx, cmap=cmapmean)
    s2.plot([img.x[1], img.x[end]], [img.y[rownum], img.y[rownum]], "--k", alpha=0.5)
    s2.set_title(central_tendency*" "*L"\log_{10}λ_x")
    s2.set_ylabel("distance y")
    s2.set_xlabel("distance x")
    # s2.set_aspect(aspect)
    s3 = subplot(233, sharex=s1, sharey=s1)
    cb2 = colorbar(im2, ax=s2)
    im3 = s3.pcolormesh(img.x, img.y, msy, cmap=cmapmean)
    s3.plot([img.x[1], img.x[end]], [img.y[rownum], img.y[rownum]], "--k", alpha=0.5)
    # s3.set_aspect(aspect)
    s3.invert_yaxis()
    s3.set_title(central_tendency*" "*L"\log_{10}λ_y")
    s3.set_ylabel("distance y")
    s3.set_xlabel("distance x")
    cb3 = colorbar(im3, ax=s3)
    s4 = subplot(234, sharex=s1)
    im4 = s4.pcolormesh(img.x, edges_r_ns[:], himage_r_ns', cmap=cmappdf)
    s4.plot(img.x, mns[rownum,:][:],"--w", alpha=0.5)
    s4.plot(img.x, CI_r_ns,"--c", alpha=0.5)
    s4.set_title("Pixel value PDF")
    s4.set_ylabel("pixel value")
    s4.set_xlabel("distance x")
    cb4 = colorbar(im4, ax=s4)
    cb4.ax.set_title("pdf")
    s5 = subplot(235, sharex=s1)
    im5 = s5.pcolormesh(img.x, edges_x_r[:], himage_x_r', cmap=cmappdf)
    cb5 = colorbar(im5, ax=s5)
    cb5.ax.set_title("pdf")
    s5.plot(img.x, msx[rownum,:][:],"--w", alpha=0.5)
    s5.plot(img.x, CI_x_r,"--c", alpha=0.5)
    s5.set_title(L"\log_{10}λ_x"*" PDF")
    s5.set_ylabel(L"\log_{10}λ_x")
    s5.set_xlabel("distance x")
    s6 = subplot(236, sharex=s1)
    im6 = s6.pcolormesh(img.x, edges_y_r[:], himage_y_r', cmap=cmappdf)
    cb6 = colorbar(im6, ax=s6)
    cb6.ax.set_title("pdf")
    s6.plot(img.x, msy[rownum,:][:],"--w", alpha=0.5)
    s6.plot(img.x, CI_y_r,"--c", alpha=0.5)
    s6.set_title(L"\log_{10}λ_y"*" PDF")
    s6.set_ylabel(L"\log_{10}λ_y")
    s6.set_xlabel("distance x")
    nicenup(f, fsize=fsize)

    f = figure(figsize=figsizecols)
    s1 = subplot(321)
    im1 = s1.pcolormesh(img.x, img.y, mns, cmap=cmapmean)
    s1.plot([img.y[colnum], img.y[colnum]], [img.y[1], img.y[end]], "--k", alpha=0.5)
    s1.set_title(central_tendency*" pixel value")
    s1.set_ylabel("distance y")
    s1.set_xlabel("distance x")
    # s1.set_aspect(aspect)
    cb1 = colorbar(im1, ax=s1)
    s2 = subplot(323, sharex=s1, sharey=s1)
    im2 = s2.pcolormesh(img.x, img.y, msx, cmap=cmapmean)
    s2.plot([img.y[colnum], img.y[colnum]], [img.y[1], img.y[end]], "--k", alpha=0.5)
    s2.set_title(central_tendency*" "*L"\log_{10}λ_x")
    s2.set_ylabel("distance y")
    s2.set_xlabel("distance x")
    # s2.set_aspect(aspect)
    s3 = subplot(325, sharex=s1, sharey=s1)
    cb2 = colorbar(im2, ax=s2)
    im3 = s3.pcolormesh(img.x, img.y, msy, cmap=cmapmean)
    s3.plot([img.y[colnum], img.y[colnum]], [img.y[1], img.y[end]], "--k", alpha=0.5)
    # s3.set_aspect(aspect)
    s3.invert_yaxis()
    s3.set_title(central_tendency*" "*L"\log_{10}λ_y")
    s3.set_ylabel("distance y")
    s3.set_xlabel("distance x")
    cb3 = colorbar(im3, ax=s3)
    s4 = subplot(322, sharey=s1)
    im4 = s4.pcolormesh(edges_c_ns[:], img.y, himage_c_ns, cmap=cmappdf)
    s4.plot(mns[:,colnum][:], img.y, "--w", alpha=0.5)
    s4.plot(CI_c_ns, img.y, "--c", alpha=0.5)
    s4.set_title("Pixel value PDF")
    s4.set_xlabel("pixel value")
    s4.set_ylabel("distance y")
    cb4 = colorbar(im4, ax=s4)
    cb4.ax.set_title("pdf")
    s5 = subplot(324, sharey=s1)
    im5 = s5.pcolormesh(edges_x_c[:], img.y, himage_x_c, cmap=cmappdf)
    cb5 = colorbar(im5, ax=s5)
    cb5.ax.set_title("pdf")
    s5.plot(msx[:,colnum][:], img.y, "--w", alpha=0.5)
    s5.plot(CI_x_c, img.y, "--c", alpha=0.5)
    s5.set_title(L"\log_{10}λ_x"*" PDF")
    s5.set_xlabel(L"\log_{10}λ_x")
    s5.set_ylabel("distance y")
    s6 = subplot(326, sharey=s1)
    im6 = s6.pcolormesh(edges_y_c[:], img.y, himage_y_c, cmap=cmappdf)
    cb6 = colorbar(im6, ax=s6)
    cb6.ax.set_title("pdf")
    s6.plot(msy[:,colnum][:], img.y, "--w", alpha=0.5)
    s6.plot(CI_y_c, img.y, "--c", alpha=0.5)
    s6.set_title(L"\log_{10}λ_y"*" PDF")
    s6.set_xlabel(L"\log_{10}λ_y")
    s6.set_ylabel("distance y")
    nicenup(f, fsize=fsize)
end

function plot_image_posterior(opt::TransD_GP.Options,
                        img::Img;
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
                        getquant = 0.5)
    @assert temperaturenum == 1
    @assert 0.0 < getquant < 1.0
    if typeof(opt) == TransD_GP.OptionsStat
        @assert opt.updatenonstat == false
        @assert opt.needλ²fromlog == false
    end
    M = assembleTat1(opt, :fstar, burninfrac=burninfrac, temperaturenum=temperaturenum)
    if domean
        m = mean(M)
    else
        a = [m' for m in M]
        m = get_image_quantile(a, getquant)
    end
    m = reshape(m, length(img.y), length(img.x))
    himage_r, edges_r, CI_r = slice_image_posterior(M, opt, img, :row, rowcolnum=rownum, nbins = nbins, qp1=qp1, qp2=qp2,
                                pdfnormalize=pdfnormalize, temperaturenum=temperaturenum)
    himage_c, edges_c, CI_c = slice_image_posterior(M, opt, img, :col, rowcolnum=colnum, nbins = nbins, qp1=qp1, qp2=qp2,
                                pdfnormalize=pdfnormalize, temperaturenum=temperaturenum)
    f = figure(figsize=figsize)
    s1 = subplot(223)
    central_tendency = domean ? ("Mean") : ("Percentile "*"$(round(Int, getquant*100))")
    im1 = s1.pcolormesh(img.x, img.y, m, cmap=cmapmean)
    s1.plot([img.x[1], img.x[end]], [img.y[rownum], img.y[rownum]], "--k", alpha=0.5)
    s1.plot([img.y[colnum], img.y[colnum]], [img.y[1], img.y[end]], "--k", alpha=0.5)
    s1.set_title(central_tendency*" pixel value")
    s1.set_ylabel("distance y")
    s1.invert_yaxis()
    # s1.set_aspect(aspect)
    s1.set_xlabel("distance x")
    cb1 = colorbar(im1, ax=s1)
    s2 = subplot(221, sharex=s1)
    im2 = s2.pcolormesh(img.x, edges_r[:], himage_r', cmap=cmappdf)
    s2.plot(img.x, m[rownum,:][:],"--w", alpha=0.5)
    s2.plot(img.x, CI_r, "--c", alpha=0.5)
    s2.set_title("Pixel value PDF")
    s2.set_ylabel("pixel value")
    s2.set_xlabel("distance x")
    cb2 = colorbar(im2, ax=s2)
    cb2.ax.set_title("pdf")
    s3 = subplot(224, sharey=s1)
    im3 = s3.pcolormesh(edges_c[:], img.y, himage_c, cmap=cmappdf)
    s3.plot(m[:,colnum][:], img.y, "--w", alpha=0.5)
    s3.plot(CI_c, img.y, "--c", alpha=0.5)
    s3.set_title("Pixel value PDF")
    s3.set_xlabel("pixel value")
    s3.set_ylabel("distance y")
    cb3 = colorbar(im3, ax=s3)
    cb3.ax.set_title("pdf")
    nicenup(f, fsize=fsize)
end

end
