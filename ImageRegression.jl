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

function getchi2forall(opt_in::TransD_GP.Options;
                        nchains          = 1,
                        figsize          = (12,6),
                        fsize            = 14,
                      )
    if nchains == 1 # then actually find out how many chains there are saved
        nchains = length(filter( x -> occursin(r"misfits.*bin", x), readdir(pwd()) )) # my terrible regex
    end
    # now look at any chain to get how many iterations
    costs_filename = "misfits_"*opt_in.fdataname
    opt_in.costs_filename    = costs_filename*"_1.bin"
    iters          = TransD_GP.history(opt_in, stat=:iter)
    niters         = length(iters)
    # then create arrays of unsorted by temperature T, k, and chi2
    Tacrosschains  = zeros(Float64, niters, nchains)
    kacrosschains  = zeros(Int, niters, nchains)
    X2by2inchains  = zeros(Float64, niters, nchains)
    # get the values into the arrays
    for ichain in 1:nchains
        opt_in.costs_filename = costs_filename*"_$ichain.bin"
        Tacrosschains[:,ichain] = TransD_GP.history(opt_in, stat=:T)
        kacrosschains[:,ichain] = TransD_GP.history(opt_in, stat=:nodes)
        X2by2inchains[:,ichain] = TransD_GP.history(opt_in, stat=:U)
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

function plot_last_target_model(img::Img, opt_in::TransD_GP.Options;
                               nchains          = 1, 
                               fsize            = 14
                               )
    if nchains == 1 # then actually find out how many chains there are saved
        nchains = length(filter( x -> occursin(r"misfits.*bin", x), readdir(pwd()) )) # my terrible regex
    end
    # now look at any chain to get how many iterations
    costs_filename = "misfits_"*opt_in.fdataname
    opt_in.costs_filename    = costs_filename*"_1.bin"
    iters          = TransD_GP.history(opt_in, stat=:iter)
    niters         = length(iters)
    # then create arrays of unsorted by temperature T
    Tacrosschains  = zeros(Float64, niters, nchains)
    # get the values into the arrays
    for ichain in 1:nchains
        opt_in.costs_filename = costs_filename*"_$ichain.bin"
        Tacrosschains[:,ichain] = TransD_GP.history(opt_in, stat=:T)
    end

    last_target_model_idx = findall(abs.(Tacrosschains[end,:] .-1.0) .< 1e-12)
    for idx in last_target_model_idx
        opt_in.fstar_filename = "models_"*opt_in.fdataname*"_$idx.bin"
        m_last = TransD_GP.history(opt_in, stat=:fstar)[end]
        figure()
        x, y = img.x, img.y
        imshow(reshape(m_last,length(y), length(x)), extent=[x[1],x[end],y[end],y[1]])
    end

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
        
