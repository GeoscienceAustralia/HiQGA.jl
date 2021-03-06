## same for all soundings
zfixed   = [-1e5]
ρfixed   = [1e12]
zstart = 0.0
extendfrac, dz = 1.06, 2.
nlayers = 40
ρbg = 10.
ntimesperdecade = 10
nfreqsperdecade = 5
## make transD options
nmin, nmax = 2, 40
K = transD_GP.GP.Mat32()
demean = true
fbounds = [-0.5 2.5]
sdpos = 0.05
sdprop = 0.05
λ, δ = [2], 0.1
save_freq = 25
Tmax = 2.50
nchainsatone = 1
## plot n random soundings and a background response
using Random
nplot = 2
nplot = min(nplot, length(sounding))
for idx in randperm(length(sounding))[1:nplot]
        aem, znall = transD_GP.SkyTEM1DInversion.makeoperator(sounding[idx],
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
                               plotfield = true)

        opt = transD_GP.SkyTEM1DInversion.make_tdgp_opt(znall = znall,
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
end
