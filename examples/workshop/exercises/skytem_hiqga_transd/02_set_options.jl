## same for all soundingss
nsamples           = 101
nchainspersounding = 4
ppn                = 5
@assert mod(ppn, nchainspersounding+1) == 0
# model discretization
zfixed             = [-1e5]
ρfixed             = [1e12]
zstart             = 0.0
showgeomplot       = true
extendfrac, dz     = 1.06, 1.25
nlayers            = 52
ρbg                = 10.
restart            = false
## make transD options
nmin, nmax         = 2, 40
K                  = transD_GP.GP.OrstUhn()
fbounds            = [-0.5 2.5] # log 10 resistivity
sdpos              = 0.05
sdprop             = 0.05
sddc               = 0.01
λ, δ               = [2], 0.1 # correlation length in depth units, tolerance for GP resistivity
save_freq          = 50
Tmax               = 2.50
## plot n random soundings and a background response
using Random
nplot = 1
nplot = min(nplot, length(soundings))
opt = transD_GP.SkyTEM1DInversion.makeoperatorandoptions(
                            rseed = 2,
                            soundings,
                            zfixed = zfixed,
                            ρfixed = ρfixed,
                            zstart = zstart,
                            extendfrac = extendfrac,
                            dz = dz,
                            ρbg = ρbg,
                            nlayers = nlayers,
                            nmin = nmin,
                            nmax = nmax,
                            K = K,
                            sdpos = sdpos,
                            sdprop = sdprop,
                            sddc = sddc,
                            fbounds = fbounds,
                            save_freq = save_freq,
                            showgeomplot = showgeomplot,
                            λ = λ,
                            δ = δ,
                            restart = restart,
                            nplot = nplot,)
