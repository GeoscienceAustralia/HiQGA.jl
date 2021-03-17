## same for all soundings
zfixed   = [-1e5]
ρfixed   = [1e12]
zstart = 0.0
extendfrac, dz = 1.06, 1.5
nlayers = 50
ρbg = 5.
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
save_freq = 50
Tmax = 2.50
nchainsatone = 1
## make nuisance options
nuisance_sdev = [
        0.0 #zTx
        0.005 #zRx
        0.005 #x_rx
        0.0 #y_rx
        0.0 #rx_roll
        0.005 #rx_pitch
        0.0 #rx_yaw
        0.0 #tx_roll
        0.0 #tx_pitch
        0.0 #tx_yaw
           ]
 nuisance_bounds = [
      0.   0.
     -5.0  5.0
     -5.0  5.0
      0.   0.
      0.   0.
     -2.0  2.0
      0.   0.
      0.   0.
      0.   0.
      0.   0.
          ]
 updatenuisances = true
## plot n random soundings and a background response
using Random
nplot = 2
nplot = min(nplot, length(sounding))
idx = 4000
# for idx in randperm(length(sounding))[1:nplot]
        aem, znall = transD_GP.TEMPEST1DInversion.makeoperator(sounding[idx],
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

        opt, optn = transD_GP.TEMPEST1DInversion.make_tdgp_opt(sounding[idx],
                                znall = znall,
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
                                nuisance_bounds = nuisance_bounds,
                                nuisance_sdev = nuisance_sdev,
                                updatenuisances = updatenuisances,
                                dispstatstoscreen = true
                                )
# end
