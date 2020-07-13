module ImageRegression
import AbstractOperator.get_misfit
using AbstractOperator
using TransD_GP, PyPlot, LinearAlgebra
using Random, Images, CommonToAll

export Img, get_image_data, calc_image_RMS, get_image_prediction_points, plot_image_data

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
                   img::Img;  s=10, fsize=14)
    f, x, y = img.f, img.x, img.y
    f1, ax1 = plt.subplots(1,2,figsize=(10,5), sharex=true, sharey=true)
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

end
