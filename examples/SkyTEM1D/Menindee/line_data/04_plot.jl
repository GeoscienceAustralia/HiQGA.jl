srcdir = dirname(dirname(dirname(dirname(pwd()))))*"/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
## plot n random soundings and a background response
using GeophysOperator, Random, PyPlot
idx = 1
## make a closure to plot posteriors
function plotposts(idx; computeforward=false, plotposterior=true, nbins=50,
                    computequants=false, nforwards=50, burninfrac=0.5, quants=[0.1,0.5,0.9])
    aem, znall = SkyTEM1DInversion.makeoperator(sounding[idx],
                               zfixed = zfixed,
                               ρfixed = ρfixed,
                               zstart = zstart,
                               extendfrac = extendfrac,
                               dz = dz,
                               ρbg = ρbg,
                               nlayers = nlayers,
                               ntimesperdecade = ntimesperdecade,
                               nfreqsperdecade = nfreqsperdecade,
                               showgeomplot = false,
                               plotfield = false)

    opt, optdummy = SkyTEM1DInversion.make_tdgp_statmode_opt(znall = znall,
                                fileprefix = sounding[idx].sounding_string,
                                nmin = nmin,
                                nmax = nmax,
                                K = K,
                                demean = demean,
                                sdpos = sdpos,
                                sdprop = sdprop,
                                fbounds = fbounds,
                                save_freq = save_freq,
                                λ = λ,
                                δ = δ,
                                )
    zall, znall, zboundaries = GeophysOperator.CommonToAll.setupz(zstart, extendfrac, dz=dz, n=nlayers)
    opt.xall[:] .= zall
    if plotposterior
        GeophysOperator.getchi2forall(opt, alpha=0.8)
        ax = gcf().axes
        ax[2].set_ylim(10,40)
        ## plot posterior
        GeophysOperator.plot_posterior(aem, opt, burninfrac=burninfrac, nbins=nbins)
        ax = gcf().axes
        ax[1].invert_xaxis()
        ## plot forward response
    end
    if (computeforward || computequants)
        M = GeophysOperator.assembleTat1(opt, :fstar, temperaturenum=1, burninfrac=burninfrac)
    end
    if computeforward
        Random.seed!(10)
        SkyTEM1DInversion.plotmodelfield!(aem, M[randperm(length(M))[1:nforwards]],
                                       dz=dz, extendfrac=extendfrac, onesigma=false, alpha=0.05)
    end
    if computequants
        m = vcat(M...)
        qs = [quantile(m[:,i],quants) for i in 1:size(m,2)]
    end
end
plotposts(idx, computeforward=true)
## now try to do this for all soundings ...
[plotposts(i, plotposterior=false, computequants=true) for i in 1:length(sounding)]
