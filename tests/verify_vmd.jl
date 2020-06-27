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
zRx = -0.02
zTx = -0.01
# Note that the receiver depth needs to be in same model layer as transmitter.
##
F = AEM_VMD_HMD.HFieldDHT(
                      nmax   = nmax,
                      zTx    = zTx,
                      rRx    = rRx,
                      freqs  = freqs,
                      zRx    = zRx)

##
AEM_VMD_HMD.getfieldFD!(F, zfixed, rho)
##
figure()
loglog(freqs,abs.(imag(F.HFD)), label="imaginary")
loglog(freqs,abs.(real(F.HFD)), label="real")
xlim(extrema(freqs))
legend()
grid()
## timing
ntimes = 1000
t = time()
for i = 1:ntimes
    AEM_VMD_HMD.getfieldFD!(F, zfixed, rho)
end
t = time() - t
@info "timing is $(t/ntimes) s"
