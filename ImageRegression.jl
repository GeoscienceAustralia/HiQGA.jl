module ImageRegression
using PyPlot, Random, Images, FileIO, DelimitedFiles, TransD_GP, MCMC_Driver


mutable struct Img

    filename          :: String
    x                 :: StepRangeLen
    y                 :: StepRangeLen
    fractrain         :: Float64
    dec               :: Int
    gausskernelwidth  :: Int 
    f                 :: Array{Float64, 2}

end

function Img(;
             filename         = "",
             dx               = 0.01,
             fractrain        = 0.02,
             dec::Int         = 2,  
             gausskernelwidth = 7)

    @assert fractrain > 0 && fractrain < 1
    @assert !(filename == "")
    @assert dec > 1
    f = get_image(filename, gausskernelwidth, dec)
    x = 0:dx:dx*size(f,2)-1
    y = 0:dx:dx*size(f,1)-1
    Img(filename, x, y, fractrain, dec, gausskernelwidth, f)
end    

function get_image(filename::String, gausskernelwidth::Int, dec)
    f = Gray.(load(filename))
    f = convert(Array{Float64, 2}, f)[1:dec:end,:1:dec:end]
    # convert image into something that looks like log resistivity
    f = -1 .+ 3*f 
    # smooth image
    imfilter(f,Kernel.gaussian(gausskernelwidth))
end    

function get_training_data(img::Img; 
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
    noisyd, δtry, noisyd[lgood], Xtrain
end

function get_all_prediction_points(img::Img)
    f, x, y = img.f, img.x, img.y
    Xall = zeros(2, length(f))
    for i in 1:length(f)
        yid, xid = Tuple(CartesianIndices(f)[i])
        Xall[:,i] = [x[xid]; y[yid]]
    end
    Xall
end

function plot_data(ftrain::Array{Float64, 1}, Xtrain::Array{Float64, 2}, 
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

function calc_simple_RMS(d::AbstractArray, img::Img, opt_in::TransD_GP.Options, 
                         opt_EM_in::MCMC_Driver.EMoptions, sd::Float64)
    
    m_true = TransD_GP.init(opt_in)
    m_true.fstar[:] = img.f
    MLnoise = opt_EM_in.MLnoise
    opt_EM_in.MLnoise = false
    @info "RMS error is" sqrt(2.0*MCMC_Driver.get_misfit(m_true, d, opt_in, opt_EM_in)/
                              sum(.!(isnan.(d))))
    opt_EM_in.MLnoise = MLnoise

end

function plot_last_target_model(img::Img, opt_in::TransD_GP.Options)
    iter_T = readdlm(opt_in.fdataname*"_temps.txt")
    last_target_model_idx = findfirst(abs.(iter_T[end,2:end] .-1.0) .< 1e-12)
    opt_in.fstar_filename = "models_"*opt_in.fdataname*"_$last_target_model_idx.bin"
    m_last = TransD_GP.history(opt_in, stat=:fstar)[end]
    figure()
    x, y = img.x, img.y
    imshow(reshape(m_last,length(y), length(x)), extent=[x[1],x[end],y[end],y[1]])
end    

function nicenup(g::PyPlot.Figure;fsize=14)
    for ax in gcf().axes
        ax.tick_params("both",labelsize=fsize)
        ax.xaxis.label.set_fontsize(fsize)
        ax.yaxis.label.set_fontsize(fsize)
        ax.title.set_fontsize(fsize)
        if typeof(ax.get_legend_handles_labels()[1]) != Array{Any,1}
            ax.legend(loc="best", fontsize=fsize)
        end
    end
    g.tight_layout()
end

end
        
