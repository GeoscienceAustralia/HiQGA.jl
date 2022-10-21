## same for all soundingss
nsamples = 200001
nchainspersounding = 5
ppn = 48
@assert mod(ppn, nchainspersounding+1) == 0
zfixed   = [-1e5]
ρfixed   = [1e12]
zstart = 0.0
showgeomplot = true
plotfield = true
extendfrac, dz = 1.06, 1.25
nlayers = 52
ρbg = 5.
ntimesperdecade = 10
nfreqsperdecade = 5
useML = false
## make transD options
nmin, nmax = 2, 40
K = transD_GP.GP.OrstUhn()
demean = false
sampledc = true
fbounds = [-1 3.0]
sdpos = 0.05
sdprop = 0.05
sddc = 0.01
λ, δ = [2], 0.1
save_freq = 100
Tmax = 2.50
nchainsatone = 1
## plot 1 random sounding and a background response
aem, opt = transD_GP.makeaemoperatorandoptions(
                            soundings[rand(1:length(soundings))];
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
                            plotfield,
                            λ = λ,
                            δ = δ
)
