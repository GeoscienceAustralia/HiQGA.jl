## same for all soundingss
nsamples = 2001
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
fbounds = [-1 3.0] # log10 resistivity
sdpos = 0.05 # fraction of fbounds
sdprop = 0.05 # fraction of nlayers
sddc = 0.01 # fraction of fbounds
λ, δ = [2], 0.1 # GP scale length in layers, and nugget in log10 
save_freq = 50
Tmax = 2.50
nchainsatone = 1
restart = false
## plot 1 random sounding and a background response
aem, opt, zall = transD_GP.makeAEMoperatorandoptions(
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
                            δ = δ, restart
)
