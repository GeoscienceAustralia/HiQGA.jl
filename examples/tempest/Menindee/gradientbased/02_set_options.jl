## for gradient descent, all model values are in log10 conductivity
σstart, σ0         = -2, -2
zfixed             = [-1e5]
ρfixed             = [1e12]
zstart             = 0.0
extendfrac, dz     = 1.06, 1.25
nlayers            = 52
ρbg                = 10
ntimesperdecade    = 10
nfreqsperdecade    = 5
regtype            = :R1
nstepsmax          = 40
ntries             = 10
lo                 = -2.5
hi                 = 0.5
λ²min              = -0.5
λ²max              = 8
β²                 = 0.0001
knownvalue         = 0.7
showgeomplot       = false
plotfield          = true
vectorsum          = true
calcjacobian       = true
## nuisance stuff
if vectorsum
    Δ = [-2.5    2.5
         -2.5    2.5]
else
    Δ = [-2.5     2.5
                -2.5     2.5
                -1       1]
end
nstepsmax      = 40 
ntriesnu       = 5
boxiters       = 3
usebox         = true
reducenuto     = 0.2
debuglevel     = 0
breaknuonknown = false
## test aem operator with one sounding    
aem, zall, = transD_GP.TEMPEST1DInversion.makeoperator(soundings[rand(1:length(soundings))];
                                    zfixed, ρfixed, zstart, extendfrac,
                                    dz, ρbg, nlayers, ntimesperdecade,
                                    nfreqsperdecade, showgeomplot,
                                    plotfield, vectorsum, calcjacobian);
