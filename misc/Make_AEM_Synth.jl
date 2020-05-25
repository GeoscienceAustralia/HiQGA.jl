module Make_AEM_Synth
using PyPlot, GP, Random, UseGA_AEM, TransD_GP

function geomprogdepth(n, dy, c) 
    dy*(1.0-c^n)/(1-c)
end

function getn(z, dy, c)
    log(1 - z/dy*(1-c))/log(c)
end    

function main()
    n, dz, extendfrac  = 52, 1.1, 1.05
    znrange            = 1.0:n
    zboundaries        = geomprogdepth.(znrange, dz, extendfrac)
    thickness          = [zboundaries[1]; diff(zboundaries)[1:end-1]]
    zall               = [zboundaries[1]/2; 0.5*(zboundaries[1:end-1] + zboundaries[2:end])]
    znall              = getn.(zall, dz, extendfrac)
    
    figure()
    plot(znall, zall)
    xlabel("depth index")
    ylabel("depth associated")
    grid()
    nicenup(gcf())    

    f, ax = plt.subplots(1, 2, figsize=(10,5))
    ax[1].stem(zboundaries[1:end-1], zboundaries[1:end-1], markerfmt="")
    ax[1].stem(zall, zall, "k--", markerfmt=" ")
    ax[1].set_xlabel("depth m")
    ax[1].set_ylabel("depth m")
    
    ax[2].stem(znrange[1:end-1], znrange[1:end-1], markerfmt="")
    ax[2].stem(znall, znall, "k--", markerfmt=" ")
    ax[2].set_ylabel("depth index")
    ax[2].set_xlabel("depth index")
    nicenup(gcf())
    
    f, ax = plt.subplots(1, 2, figsize=(10,5), sharey=true)
    ax[1].stem(zall[1:end-1],thickness, "k--", markerfmt=" ")
    ax[1].set_ylabel("thickness m")
    ax[1].yaxis.grid(which="major")
    ax[2].stem(znall[1:end-1],thickness, "k--", markerfmt=" ")
    ax[2].set_xlabel("depth index")
    ax[2].yaxis.grid(which="major")
    nicenup(f)

    λ, δ = [2.0], 0.01
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
