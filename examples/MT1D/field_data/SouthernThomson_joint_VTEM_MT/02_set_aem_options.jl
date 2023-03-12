showgeomplot = false
plotfield = true
ntimesperdecade = 10
nfreqsperdecade = 5
useML = false
ρbg = 10.
## make transD options
nmin, nmax = 2, 50
K = transD_GP.GP.OrstUhn()
fbounds = [-1 3.0]
sdpos = 0.05
sdprop = 0.05
sddc = 0.01
λ, δ = [2], 0.1
save_freq = 25
Tmax = 2.50
nchainsatone = 1
## plot 1 random sounding and a background response
aem, opt, zall = transD_GP.makeAEMoperatorandoptions(
                            soundings[1];
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
                            sdpos = sdpos,
                            sdprop = sdprop,
                            sddc = sddc,
                            fbounds = fbounds,
                            save_freq = save_freq,
                            showgeomplot = showgeomplot,
                            plotfield,
                            λ = λ,
                            δ = δ
)
opt.dispstatstoscreen = false