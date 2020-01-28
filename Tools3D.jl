module Tools3D

using PyPlot, Random, TransD_GP, GP, Statistics

function geomprogdepth(n, dy, c)
    dy*(1.0-c^n)/(1-c)
end

function getn(z, dy, c)
    log(1 - z/dy*(1-c))/log(c)
end

function nicenup(g::PyPlot.Figure;fsize=16)
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

function meshgrid(xs, ys)
    xx, yy = [xs[i] for i in 1:length(xs), j in 1:length(ys)], [ys[j] for i in 1:length(xs), j in 1:length(ys)]
    return xx, yy
end

function makezn(;n          = 20,
                 dz         = 1.6,
                 extendfrac = 1.245)
    @assert typeof(n) == Int
    @assert dz > 0.0
    @assert extendfrac > 1.0
    znrange            = 1.0:n
    zboundaries        = geomprogdepth.(znrange, dz, extendfrac)
    thickness          = [zboundaries[1]; diff(zboundaries)[1:end-1]]
    zall               = [zboundaries[1]/2; 0.5*(zboundaries[1:end-1] + zboundaries[2:end])]
    znall              = getn.(zall, dz, extendfrac)

    figure()
    plot(znall, zall)
    xlabel("depth index")
    ylabel("depth associated km")
    grid()
    nicenup(gcf())

    f, ax = plt.subplots(1, 2, figsize=(10,5))
    ax[1].stem(zboundaries[1:end-1], zboundaries[1:end-1], markerfmt="")
    ax[1].stem(zall, zall, "k--", markerfmt=" ")
    ax[1].set_xlabel("depth km")
    ax[1].set_ylabel("depth km")

    ax[2].stem(znrange[1:end-1], znrange[1:end-1], markerfmt="")
    ax[2].stem(znall, znall, "k--", markerfmt=" ")
    ax[2].set_ylabel("depth index")
    ax[2].set_xlabel("depth index")
    nicenup(gcf())

    f, ax = plt.subplots(1, 2, figsize=(10,5), sharey=true)
    ax[1].stem(zall[1:end-1],thickness, "k--", markerfmt=" ")
    ax[1].set_xlabel("depth km")
    ax[1].set_ylabel("thickness km")
    ax[1].yaxis.grid(which="major")
    ax[2].stem(znall[1:end-1],thickness, "k--", markerfmt=" ")
    ax[2].set_xlabel("depth index")
    ax[2].yaxis.grid(which="major")
    nicenup(f)
    return zall, znall
end

function makeopt(;nmin      = 2,
                  nmax      = 200,
                  λ         = [150.0, 400.0, 1.0],
                  δ         = 0.1,
                  fbounds   = [0.5 4],
                  demean    = true,
                  sdev_prop = 0.1,
                  sdev_pos  = [8;8.0;1.4],
                  pnorm     = 2.0,
                  λx        = 1600.0,
                  λy        = 1600.0,
                  dxfrac    = 0.025,
                  dyfrac    = 0.025,
                  rseed     = 4,
                  timebirth = true,
                  znall     = nothing
                )
    @assert znall != nothing
    dx = dxfrac*λx
    dy = dyfrac*λy
    x = 0:(dx):λx-dx
    y = 0:(dy):λy-dy
    xall = zeros(3,length(x)*length(y)*length(znall))
    for i in 1:size(xall,2)
        xid, yid, zid = Tuple(CartesianIndices((length(x),length(y),length(znall)))[i])
        xall[:,i] = [x[xid]; y[yid]; znall[zid]]
    end
    xbounds = zeros(Float64,size(xall,1),2)
    for dim in 1:size(xall, 1)
        xbounds[dim,:] = [minimum(xall[dim,:]), maximum(xall[dim,:])]
    end
    ## Initialize a model using these options
    Random.seed!(rseed)
    opt = TransD_GP.Options(nmin = nmin,
                            nmax = nmax,
                            xbounds = xbounds,
                            fbounds = fbounds,
                            xall = xall,
                            λ = λ,
                            δ = δ,
                            demean = demean,
                            sdev_prop = sdev_prop,
                            sdev_pos = sdev_pos,
                            pnorm = pnorm,
                            quasimultid = false
                            )
    return opt, x, y
end

function slicemodel(m::TransD_GP.Model,
                    opt::TransD_GP.Options;
                    slicesx    = nothing,
                    slicesy    = nothing,
                    slicesz    = nothing,
                    dz         = nothing,
                    extendfrac = nothing,
                    )
    @assert !any((slicesx, slicesy, slicesz,
                  dz, extendfrac) .== nothing)
    @assert !(slicesz != [] && (slicesx != [] || slicesy != []))
    @assert !(slicesx != [] && (slicesz != [] || slicesy != []))
    @assert !(slicesx == [] && slicesy == [] && slicesz ==[])
    slicemodel(m.fstar, m.n, m.xtrain, m.ftrain, opt::TransD_GP.Options; slicesx=slicesx, slicesy=slicesy, slicesz=slicesz,
                dz=dz, extendfrac=extendfrac)
end

function slicemodel(fstar::Array{Float64, 1}, npoints::Int,
                    xtrain::Array{Float64, 2}, ftrain::Array{Float64, 1},
                    opt::TransD_GP.Options;
                    slicesx    = nothing,
                    slicesy    = nothing,
                    slicesz    = nothing,
                    dz         = nothing,
                    extendfrac = nothing,
                    )
    x = unique(opt.xall[1,:])
    y = unique(opt.xall[2,:])
    z = geomprogdepth.(unique(opt.xall[3,:]), dz, extendfrac)
    @assert !any((slicesx, slicesy, slicesz,
                  dz, extendfrac) .== nothing)
    @assert !(slicesz != [] && (slicesx != [] || slicesy != []))
    @assert !(slicesx != [] && (slicesz != [] || slicesy != []))
    @assert !(slicesx == [] && slicesy == [] && slicesz ==[])
    v = reshape(fstar,length(x), length(y), length(z))
    f = figure(figsize=(10,10))
    if slicesz!=[]
        xx,yy = meshgrid(x,y)
        l = slicesz
        for i in l
            c=(v[:,:,i].-minimum(v))./(maximum(v)-minimum(v))
            plot_surface(xx, yy, z[i]*ones(size(yy)), facecolors=plt.cm.jet_r(c), shade=false)
        end
    elseif slicesx!=[]
        yy,zz = meshgrid(y,z)
        l = slicesx
        for i in l
            c=(v[i,:,:].-minimum(v))./(maximum(v)-minimum(v))
            plot_surface(x[i]*ones(size(yy)), yy, zz, facecolors=plt.cm.jet_r(c), shade=false)
        end
    elseif slicesy!=[]
        xx,zz = meshgrid(x,z)
        l = slicesy
        for i in l
            c=(v[:,i,:].-minimum(v))./(maximum(v)-minimum(v))
            plot_surface(xx, y[i]*ones(size(zz)), zz, facecolors=plt.cm.jet_r(c), shade=false)
        end
    end
    if npoints>0
        scatter3D(xtrain[1,1:npoints], xtrain[2,1:npoints], geomprogdepth.(xtrain[3,1:npoints], dz, extendfrac),
                c=ftrain[1:npoints], vmin=minimum(v), vmax=maximum(v), s=50, cmap="jet_r")
    end
    xlim(extrema(x))
    ylim(extrema(y))
    zlim(extrema(z))
    xlabel("x km")
    ylabel("y km")
    zlabel("depth km")
    cbar = f.colorbar(plt.cm.ScalarMappable(cmap="jet_r",norm=matplotlib.colors.Normalize(vmin=minimum(v), vmax=maximum(v))),ax=gca())
    cbar.set_label(L"\log_{10} \rho")
    gca().zaxis.label.set_fontsize(14); nicenup(gcf())
    gca().invert_zaxis()
end

function makecubemodel(opt::TransD_GP.Options;
                        dz         = nothing,
                        extendfrac = nothing,
                        z1         = 10.0,
                        z2         = 400.0,
                        yleft      = 500.0,
                        yright     = 1100.0,
                        xup        = 980.0,
                        xdown      = 580.0,
                        zup        = 70.0,
                        zdown      = 250.0,
                        ρ0         = 2,
                        ρ1         = 4,
                        ρ2         = 0.5,
                        ρanom      = 0.8
                        )
    @assert extendfrac != nothing
    @assert dz != nothing
    @assert z1<zup<zdown<z2
    @assert yleft<yright
    @assert xdown<xup

    ρ = zeros(size(opt.xall, 2))
    zmin, zmax = geomprogdepth.(extrema(opt.xall[3,:]), dz, extendfrac)
    for i = 1:size(opt.xall, 2)
        x, y, zn = opt.xall[:,i]
        z = geomprogdepth(zn, dz, extendfrac)
        if z<z1
            ρ[i] = ρ0
        elseif (zup<z<zdown) & (yleft<y<yright) & (xdown<x<xup)
            ρ[i] = ρanom*(ρ1 + (z-zmin)/(zmax-z1)*(ρ2-ρ1)) #ρanom
        elseif z>z2
            ρ[i] = ρ2
        else
            ρ[i] = ρ1 + (z-zmin)/(zmax-z1)*(ρ2-ρ1)
        end
    end
    opt.fbounds = [0.98*minimum(ρ) 1.02*maximum(ρ)]
    return ρ
end

function get_training_data(fstar::Array{Float64, 1},
                           opt::TransD_GP.Options;
                           dz = nothing,
                           extendfrac = nothing,  
                           sdmaxfrac = 0.05,
                           rseed     = 12,
                           zbreak    = 12321.0,
                           takeevery = 4,
                           fractrain = 0.05 
                          )
    @assert extendfrac != nothing
    @assert dz != nothing
    @assert zbreak != 12321.0
    @assert sdmaxfrac > 0 && sdmaxfrac < 1
    Random.seed!(rseed)
    x, y =  unique(opt.xall[1,:]), unique(opt.xall[2,:])
    z = geomprogdepth.(unique(opt.xall[3,:]), dz, extendfrac)
    @assert minimum(z) < zbreak < maximum(z)
    noisyd = NaN .+ zeros(Float64, size(fstar))
    ntrain = round(Int, fractrain*length(fstar))
    linidx = randperm(length(fstar))[1:ntrain]
    δtry = sdmaxfrac*maximum(abs.(fstar))
    for l in linidx
        row, col, zid = Tuple(CartesianIndices((length(x),length(y),length(z)))[l])
        if z[zid] > zbreak
            noisyd[l] = fstar[l] + δtry*randn()
        else
            if rem(row, takeevery) == 0 && rem(col, takeevery) == 0 && rem(zid, takeevery) == 0
                noisyd[l] = fstar[l] + δtry*randn()
            end
        end
    end
    f = figure(figsize=(10,10))
    scatter3D(opt.xall[1,:], opt.xall[2,:], geomprogdepth.(opt.xall[3,:], dz, extendfrac), 
                          c=noisyd, vmin=minimum(fstar), vmax=maximum(fstar), cmap="jet_r")
    xlim(extrema(x))
    ylim(extrema(y))
    zlim(extrema(z))
    xlabel("x km")
    ylabel("y km")
    zlabel("depth km")
    cbar = f.colorbar(plt.cm.ScalarMappable(cmap="jet_r",norm=matplotlib.colors.Normalize(vmin=minimum(fstar), vmax=maximum(fstar))),ax=gca())
    cbar.set_label(L"\log_{10} \rho")
    gca().zaxis.label.set_fontsize(14); nicenup(gcf())
    gca().invert_zaxis()

    noisyd
end

function plot_last_target_model(m::TransD_GP.Model,
                                opt_in::TransD_GP.Options;
                                slicesx    = nothing,
                                slicesy    = nothing,
                                slicesz    = nothing,
                                dz         = nothing,
                                extendfrac = nothing,
                                nchains    = 1, 
                                fsize      = 14
                               )
    @assert !any((slicesx, slicesy, slicesz,
                  dz, extendfrac) .== nothing)
    @assert !(slicesz != [] && (slicesx != [] || slicesy != []))
    @assert !(slicesx != [] && (slicesz != [] || slicesy != []))
    @assert !(slicesx == [] && slicesy == [] && slicesz ==[])

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
        slicemodel(m, opt, slicesx=slicesx, slicesy=slicesy, slicesz=slicesz, dz=dz, extendfrac=extendfrac)        
    end

end

end
