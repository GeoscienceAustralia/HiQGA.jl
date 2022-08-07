## make gradientinv options same for all soundings
σstart, σ0         = -2, -2
zfixed             = [-1e5]
ρfixed             = [1e12]
zstart             = 0.0
extendfrac, dz     = 1.06, 1.25
nlayers            = 52
ρbg                = 10
ntimesperdecade    = 10
nfreqsperdecade    = 5
modelprimary       = false
regtype            = :R1
nstepsmax          = 40
ntries             = 6
target             = nothing
lo                 = -2.5
hi                 = 0.5
λ²min              = -0.5
λ²max              = 8
λ²frac             = 4
ntestdivsλ²        = 50
αmin               = -4 
αmax               = 0 
αfrac              = 4
β²                 = 0.1
ntestdivsα         = 32
regularizeupdate   = false
knownvalue         = 0.7
firstvalue         = :last
κ                  = transD_GP.GP.Mat52()
breakonknown       = true
dobo               = false
## plot n random soundings and a background response
using Random
nplot = 2
nplot = min(nplot, length(soundings))
opt = transD_GP.SkyTEM1DInversion.makeoperatorandoptions(
                            soundings,
                            rseed = 2,
                            zfixed = zfixed,
                            ρfixed = ρfixed,
                            zstart = zstart,
                            extendfrac = extendfrac,
                            dz = dz,
                            ρbg = ρbg,
                            nlayers = nlayers,
                            ntimesperdecade = ntimesperdecade,
                            nfreqsperdecade = nfreqsperdecade)

