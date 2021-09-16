using PyPlot, transD_GP, Random
## frequencies
nFreqsPerDecade     = 7
freqLowLimit        = 1e-1
freqHighLimit       = 1e5
freqs = 10 .^(log10(freqLowLimit):1/nFreqsPerDecade:log10(freqHighLimit))
## model
zfixed   = [-1e5,   0,   (20:20:200)...]
rho      = [1e12,   100, 100*ones(10)...]
rho[5] = 1
##  geometry
rRx = 100.
zRx = -37.
zTx = -35.
nkᵣeval = 200
## VMD
modelprimary = true
Fvmd = transD_GP.AEM_VMD_HMD.HFieldDHT(
                      calcjacobian = true,
                      nmax = length(zfixed),
                      zTx    = zTx,
                      rRx    = rRx,
                      freqs  = freqs,
                      zRx    = zRx,
                      nkᵣeval = nkᵣeval,
                      modelprimary = modelprimary)
transD_GP.AEM_VMD_HMD.getfieldFD!(Fvmd, zfixed, rho)