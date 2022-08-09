# model discretization
zfixed             = [-1e5]
ρfixed             = [1e12]
zstart             = 0.0
extendfrac, dz     = 1.06, 1.25
nlayers            = 52
# make gradientinv options same for all soundings
σstart, σ0         = -2, -2
ρbg                = 10 # test linear resistivity for forward model before inversion
regtype            = :R1
nstepsmax          = 10 # total iterations of Occam (outer)
ntries             = 6 # tries per outer iteration
target             = nothing
lo                 = -2.5 # log10 conductivity min
hi                 = 0.5  # log10 conductivity max
λ²min              = -0.5
λ²max              = 8
β²                 = 0.1
knownvalue         = 0.7 # exit inner loop at this much misfit reduction
breakonknown       = true # exit inner loop at knownvalue
showgeomplot       = true # show model discretization
## plot n random soundings and a background response
using Random
nplot = 1
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
                            nplot = nplot,
                            showgeomplot = showgeomplot);

