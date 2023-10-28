using PyPlot, Revise, HiQGA.transD_GP, Random
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
doconvramp = false
## VMD
modelprimary = true
Fvmd = transD_GP.AEM_VMD_HMD.HFieldDHT(;
                      nmax   = nmax,
                      zTx    = zTx,
                      rRx    = rRx,
                      doconvramp,
                      freqs  = freqs,
                      zRx    = zRx,
                      nkᵣeval = nkᵣeval,
                      modelprimary = modelprimary)
transD_GP.AEM_VMD_HMD.getfieldFD!(Fvmd, zfixed, rho)
## now use a tiny loop radius - VMD approximation is worse as radius gets larger
rTx = 0.001
Floop = transD_GP.AEM_VMD_HMD.HFieldDHT(;
                      nmax   = nmax,
                      zTx    = zTx,
                      rRx    = rRx,
                      rTx    = rTx,
                      freqs  = freqs,
                      doconvramp,
                      zRx    = zRx,
                      nkᵣeval = nkᵣeval,
                      modelprimary = modelprimary)
transD_GP.AEM_VMD_HMD.getfieldFD!(Floop, zfixed, rho)
## Compare H vertical
f, ax = plt.subplots(2,1, sharex=true)
ax[1].loglog(freqs,abs.(imag(Fvmd.HFD_z)), label="imaginary, VMD")
ax[1].loglog(freqs,abs.(real(Fvmd.HFD_z)), label="real, VMD")
ax[1].loglog(freqs,abs.(imag(Floop.HFD_z)), "*", label="imaginary, loop")
ax[1].loglog(freqs,abs.(real(Floop.HFD_z)), "*", label="real, loop")
ax[1].legend()
ax[1].set_title("Compare with W&H Fig 4.2")
## Compare H radial
getradialH = true
Fvmd = transD_GP.AEM_VMD_HMD.HFieldDHT(;
                      nmax   = nmax,
                      zTx    = 0,
                      rRx    = rRx,
                      freqs  = freqs,
                      zRx    = 0,
                      doconvramp,
                      nkᵣeval = nkᵣeval,
                      modelprimary = false, # says W&H for within plane of observation !!
                      getradialH = getradialH)
transD_GP.AEM_VMD_HMD.getfieldFD!(Fvmd, zfixed, rho)
ax[2].loglog(freqs,abs.(imag(Fvmd.HFD_r)), label="imaginary, VMD")
ax[2].loglog(freqs,abs.(real(Fvmd.HFD_r)), label="real, VMD")
ax[2].set_title("Compare with W&H Fig 4.3")
ax[2].set_xlim(extrema(freqs))
ax[2].set_ylim(1e-12,1e-6)
ax[1].grid()
ax[2].grid()
ax[2].set_xlabel("Hz")
legend()
## now test loop field, receiver at center
zTx = -0.01
zRx = -0.02
rRx = 0.001
rTx = 50.0
Floopin = transD_GP.AEM_VMD_HMD.HFieldDHT(;
                      nmax   = nmax,
                      zTx    = zTx,
                      rRx    = rRx,
                      rTx    = rTx,
                      freqs  = freqs,
                      zRx    = zRx,
                      doconvramp,
                      nkᵣeval = nkᵣeval,
                      modelprimary = modelprimary)
transD_GP.AEM_VMD_HMD.getfieldFD!(Floopin, zfixed, rho)
##
figure()
loglog(freqs,pi*rTx^2*abs.(imag(Floopin.HFD_z)), label="imaginary, rx in = $(Floopin.rxwithinloop)")
loglog(freqs,pi*rTx^2*abs.(real(Floopin.HFD_z)), label="real, rx in = $(Floopin.rxwithinloop)")
xlim(extrema(freqs))
xlabel("Hz")
legend()
title("Compare with W&H Fig 4.7")
grid()
