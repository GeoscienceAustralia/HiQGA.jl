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
zRx = -0.011
zTx = -0.01
nkᵣeval = 200
## VMD
modelprimary = true
Fvmd = AEM_VMD_HMD.HFieldDHT(
                      nmax   = nmax,
                      zTx    = zTx,
                      rRx    = rRx,
                      freqs  = freqs,
                      zRx    = zRx,
                      nkᵣeval = nkᵣeval,
                      modelprimary = modelprimary)
AEM_VMD_HMD.getfieldFD!(Fvmd, zfixed, rho)
## now use a tiny loop radius - VMD approximation is worse as radius gets larger
rTx = 0.001
Floop = AEM_VMD_HMD.HFieldDHT(
                      nmax   = nmax,
                      zTx    = zTx,
                      rRx    = rRx,
                      rTx    = rTx,
                      freqs  = freqs,
                      zRx    = zRx,
                      nkᵣeval = nkᵣeval,
                      modelprimary = modelprimary)
AEM_VMD_HMD.getfieldFD!(Floop, zfixed, rho)
## Compare H vertical
f, ax = plt.subplots(2,1, sharex=true)
ax[1].loglog(freqs,abs.(imag(Fvmd.HFD_z)), label="imaginary, VMD")
ax[1].loglog(freqs,abs.(real(Fvmd.HFD_z)), label="real, VMD")
ax[1].loglog(freqs,abs.(imag(Floop.HFD_z)), "*", label="imaginary, loop")
ax[1].loglog(freqs,abs.(real(Floop.HFD_z)), "*", label="real, loop")
ax[1].legend()
ax[1].set_title("Compare with W&H Fig 4.2")
## Compare H radial
modelprimary = false # says W&H for within plane of observation !!
getradialH = true
Fvmd = AEM_VMD_HMD.HFieldDHT(
                      nmax   = nmax,
                      zTx    = zTx,
                      rRx    = rRx,
                      freqs  = freqs,
                      zRx    = zRx,
                      nkᵣeval = nkᵣeval,
                      modelprimary = modelprimary,
                      getradialH = getradialH)
AEM_VMD_HMD.getfieldFD!(Fvmd, zfixed, rho)
ax[2].loglog(freqs,abs.(imag(Fvmd.HFD_r)), label="imaginary, VMD")
ax[2].loglog(freqs,abs.(real(Fvmd.HFD_r)), label="real, VMD")
ax[2].set_title("Compare with W&H Fig 4.3")
ax[2].set_xlim(extrema(freqs))
ax[1].grid()
ax[2].grid()
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
loglog(freqs,pi*rTx^2*abs.(imag(Floopin.HFD_z)), label="imaginary, rx in = $(Floopin.rxwithinloop)")
loglog(freqs,pi*rTx^2*abs.(real(Floopin.HFD_z)), label="real, rx in = $(Floopin.rxwithinloop)")
xlim(extrema(freqs))
legend()
title("Compare with W&H Fig 4.7")
grid()
