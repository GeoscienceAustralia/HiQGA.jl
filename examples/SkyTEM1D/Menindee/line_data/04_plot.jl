srcdir = dirname(dirname(dirname(dirname(pwd()))))*"/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
## plot n random soundings and a background response
using GeophysOperator, Random, PyPlot
idx = 200
##
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
## get stats
GeophysOperator.getchi2forall(opt, alpha=0.8)
ax = gcf().axes
ax[2].set_ylim(10,40)
## plot posterior
zall, znall, zboundaries = GeophysOperator.CommonToAll.setupz(zstart, extendfrac, dz=dz, n=nlayers)
opt.xall[:] .= zall
GeophysOperator.plot_posterior(aem, opt)
## plot forward response
M = GeophysOperator.assembleTat1(opt, :fstar, temperaturenum=1)
Random.seed!(10)
SkyTEM1DInversion.plotmodelfield!(aem, M[randperm(length(M))[1:50]],
                                       dz=dz, extendfrac=extendfrac, onesigma=false, alpha=0.05)
