srcdir = dirname(pwd())*"/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
using PyPlot, Revise, AEM_VMD_HMD, Random
## set up
nfreqsperdecade     = 6
ntimesperdecade     = 10
modelprimary        = true
provideddt          = true
doconvramp          = false
nkᵣeval             = 200
times               = 10 .^LinRange(-5,-1, 50)
lowpassfcs          = [1e6]
## model
zfixed   = [-1e5,   0,    ]
rho      = [1e12,   100,  ]
nmax = 200
##  geometry
rRx = 100.
zRx = 0.
zTx = 0.
modelprimary = true
##
F = AEM_VMD_HMD.HFieldDHT(
                      zTx    = zTx,
                      rRx    = rRx,
                      times = times,
                      zRx    = zRx,
                      nkᵣeval = nkᵣeval,
                      modelprimary = modelprimary,
                      provideddt = provideddt,
                      lowpassfcs = lowpassfcs)
AEM_VMD_HMD.getfieldTD!(F, zfixed, rho)
## plot FD
figure()
loglog(F.freqs, abs.(real.(F.HFD_z)))
loglog(F.freqs, abs.(imag.(F.HFD_z)))
xlim(1e-1,1e5)
ylim(1e-12, 1e-6)
## plot time domain
figure()
loglog(F.interptimes*1e3, abs.(F.HTD_z_interp))
xlim(1e-5,1e3)
ylim(1e-11,1e-1)
grid(true, which="both")
