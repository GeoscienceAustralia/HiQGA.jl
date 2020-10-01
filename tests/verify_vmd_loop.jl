srcdir = dirname(pwd())*"/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
using PyPlot, Revise, AEM_VMD_HMD, Random
## frequencies
nFreqsPerDecade     = 7
freqLowLimit        = 1e-1
freqHighLimit       = 1e5
freqs = 10 .^(log10(freqLowLimit):1/nFreqsPerDecade:log10(freqHighLimit))
## model
zfixed   = [-1e5,   0,    ]
rho      = [1e12,   100,  ]
nmax = 200
##  geometry
rRx = 100.
zRx = -0.01
zTx = -0.02
modelprimary = true
nkᵣeval = 150
# Note that the receiver depth needs to be in same model layer as transmitter.
##
Fvmd = AEM_VMD_HMD.HFieldDHT(
                      nmax   = nmax,
                      zTx    = zTx,
                      rRx    = rRx,
                      freqs  = freqs,
                      zRx    = zRx,
                      nkᵣeval = nkᵣeval,
                      modelprimary = modelprimary)

## now use a tiny loop radius - VMD approximation is worse as radius gets larger
##
rTx = 0.01
Floop = AEM_VMD_HMD.HFieldDHT(
                      nmax   = nmax,
                      zTx    = zTx,
                      rRx    = rRx,
                      rTx    = rTx,
                      freqs  = freqs,
                      zRx    = zRx,
                      nkᵣeval = nkᵣeval,
                      modelprimary = modelprimary)
##
AEM_VMD_HMD.getfieldFD!(Fvmd, zfixed, rho)
AEM_VMD_HMD.getfieldFD!(Floop, zfixed, rho)
##
figure()
loglog(freqs,abs.(imag(Fvmd.HFD)), label="imaginary, VMD")
loglog(freqs,abs.(real(Fvmd.HFD)), label="real, VMD")
loglog(freqs,abs.(imag(Floop.HFD)), "*", label="imaginary, loop")
loglog(freqs,abs.(real(Floop.HFD)), "*", label="real, loop")
xlim(extrema(freqs))
legend()
title("Compare with W&H Fig 4.2")
grid()
## now test loop field, receiver at center
zTx = -0.01
zRx = -0.02
rRx = 0.001
rTx = 50.0
Floopin = AEM_VMD_HMD.HFieldDHT(
                      nmax   = nmax,
                      zTx    = zTx,
                      rRx    = rRx,
                      rTx    = rTx,
                      freqs  = freqs,
                      zRx    = zRx,
                      nkᵣeval = nkᵣeval,
                      modelprimary = modelprimary)
AEM_VMD_HMD.getfieldFD!(Floopin, zfixed, rho)
##
figure()
loglog(freqs,pi*rTx^2*abs.(imag(Floopin.HFD)), label="imaginary, rx in = $(Floopin.rxwithinloop)")
loglog(freqs,pi*rTx^2*abs.(real(Floopin.HFD)), label="real, rx in = $(Floopin.rxwithinloop)")
xlim(extrema(freqs))
legend()
title("Compare with W&H Fig 4.7")
grid()
