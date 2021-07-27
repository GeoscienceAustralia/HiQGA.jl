## same for all soundingss
nsamples = 100001
nchainspersounding = 5
ppn = 48
@assert mod(ppn, nchainspersounding+1) == 0
zfixed   = [-1e5]
ρfixed   = [1e12]
zstart = 0.0
showgeomplot = false
extendfrac, dz = 1.06, 1.5
nlayers = 50
ρbg = 10.
ntimesperdecade = 10
nfreqsperdecade = 5
useML = false
## make transD options
nmin, nmax = 2, 40
K = transD_GP.GP.OrstUhn()
demean = false
sampledc = true
fbounds = [-1 3]
sdpos = 0.05
sdprop = 0.05
sddc = 0.008
λ, δ = [2], 0.1
save_freq = 100
Tmax = 2.50
nchainsatone = 1
## plot n random soundings and a background response
using Random
nplot = 2
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
                            ntimesperdecade = ntimesperdecade,
                            nfreqsperdecade = nfreqsperdecade,
                            useML = useML,
                            nmin = nmin,
                            nmax = nmax,
                            K = K,
                            demean = demean,
                            sdpos = sdpos,
                            sdprop = sdprop,
                            sddc = sddc,
                            sampledc = sampledc,
                            fbounds = fbounds,
                            save_freq = save_freq,
                            showgeomplot = showgeomplot,
                            λ = λ,
                            δ = δ
)
