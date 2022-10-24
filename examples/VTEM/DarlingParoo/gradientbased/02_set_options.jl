## make gradientinv options same for all soundings
σstart, σ0         = -2, -2
zfixed             = [-1e5]
ρfixed             = [1e12]
zstart             = 0.0
extendfrac         = 1.06
dz                 = 1.5
ρbg                = 10
nlayers            = 50
ntimesperdecade    = 10
nfreqsperdecade    = 5
regtype            = :R1
nstepsmax          = 40
ntries             = 6
lo                 = -3.
hi                 = 1.
λ²min              = -0.5
λ²max              = 8
β²                 = 0.01
knownvalue         = 0.7
breakonknown       = true
calcjacobian       = true
## plot a random sounding and a background response
aem, zall, = transD_GP.VTEM1DInversion.makeoperator(
    soundings[rand(1:length(soundings))];
    zfixed, ρfixed, zstart, extendfrac, calcjacobian,
    dz, ρbg, nlayers, ntimesperdecade, nfreqsperdecade, plotfield=true)