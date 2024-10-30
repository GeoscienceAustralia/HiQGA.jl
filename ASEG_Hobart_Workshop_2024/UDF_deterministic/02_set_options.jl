# # Occam's inversion
## make gradient inversion options same for all soundings
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
modelprimary       = false
regtype            = :R1
nstepsmax          = 40
ntries             = 10
target             = nothing
lo                 = -3. # log 10 S/m
hi                 = 1.  # log 10 S/m
λ²min              = -0.5
λ²max              = 8
λ²frac             = 4
β²                 = 0.1
knownvalue         = 0.7
breakonknown       = true
zipsaveprefix      = basename(pwd())*"_β²_$(β²)_$(regtype)_bg_$(round(10. ^σ0, sigdigits=4))_Spm"
## plot a random sounding and a background response
aem, zall, = transD_GP.SkyTEM1DInversion.makeoperator(
    soundings[rand(1:length(soundings))];
    zfixed, ρfixed, zstart, extendfrac, calcjacobian=true,
    dz, ρbg, nlayers, ntimesperdecade, nfreqsperdecade, plotfield=true)
    


