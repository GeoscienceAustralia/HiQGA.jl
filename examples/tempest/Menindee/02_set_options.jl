## same for all soundingss
nsamples = 100001
nchainspersounding = 7
ppn = 48
@assert mod(ppn, nchainspersounding+1) == 0
zfixed   = [-1e5]
ρfixed   = [1e12]
zstart = 0.0
extendfrac, dz = 1.06, 1.25
nlayers = 52
ρbg = 5.
ntimesperdecade = 10
nfreqsperdecade = 5
## make transD options
nmin, nmax = 2, 40
K = transD_GP.GP.OrstUhn()
demean = false
sampledc = true
fbounds = [-1.0 3.0]
sddc = 0.01 # fraction
sdpos = 0.05 # fraction
sdprop = 0.05 # fraction
λ, δ = [2.], 0.1
save_freq = 50
Tmax = 2.50
nchainsatone = 1
vectorsum = true
useML = false
## make nuisance options
# sdev are fractions of width
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
     -3.0  3.0
     -3.0  3.0
      0.   0.
      0.   0.
     -1.5  1.5
      0.   0.
      0.   0.
      0.   0.
      0.   0.
          ]
updatenuisances = true
## see if options set can make a reasonable looking forward
transD_GP.TEMPEST1DInversion.makeoperatorandoptions(
        soundings;
        zfixed             = zfixed,
        ρfixed             = ρfixed,
        zstart             = zstart,
        extendfrac         = extendfrac,
        dz                 = dz,
        ρbg                = ρbg,
        nlayers            = nlayers,
        ntimesperdecade    = ntimesperdecade,
        nfreqsperdecade    = nfreqsperdecade,
        nmin               = nmin,
        nmax               = nmax,
        K                  = K,
        demean             = demean,
        sampledc           = sampledc,
        sddc               = sddc,
        sdpos              = sdpos,
        sdprop             = sdprop,
        fbounds            = fbounds,
        save_freq          = save_freq,
        λ                  = λ,
        δ                  = δ,
        nuisance_bounds    = nuisance_bounds,
        nuisance_sdev      = nuisance_sdev,
        updatenuisances    = updatenuisances,
        useML              = useML,
        vectorsum          = vectorsum)