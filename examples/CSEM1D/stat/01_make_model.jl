using PyPlot, DelimitedFiles, Random, Revise,
      transD_GP
##
Random.seed!(23)
rMin = 500 #m
rMax = 5200  #m
nRx  = 32
freqs   = [0.25, 0.75, 1.75, 3.25] #Hz
rRx    = collect(LinRange(rMin,rMax,nRx))  # Ranges to receivers   (m)
zRx    = 1000-eps(1.)
zTx    = 975.
RxAzim = 0.
TxDip  = 0.
# Note that the receiver depth needs to be in same model layer as transmitter.
zfixed = [-1e6,    0      ]
ρfixed = [1e13,    0.3    ]
zstart = zRx + eps(1.)
extendfrac, dz = 1.008, 12.
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=126)
z, ρ, nfixed = transD_GP.makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
##
F = transD_GP.CSEM1DInversion.CSEM1DEr.RadialErLagged(zTx    = [zTx],
                      rRx    = rRx,
                      freqs  = freqs,
                      zRx    = [zRx],
                      RxAzim = RxAzim,
                      TxDip  = TxDip)
##
ρ[(z.>=zstart) .& (z.<1400)] .= 10^0.0
ρ[(z.>=1400) .&(z.<1700)] .= 10^0.3
ρ[(z.>=1700) .&(z.<1800)] .= 10^0.5
ρ[(z.>=1800) .&(z.<1900)] .= 10^0.2
ρ[(z.>=1900) .&(z.<2000)] .= 10^1.4 # 1.4 or 0.5
ρ[(z.>=2000) .&(z.<2400)] .= 10^0.0
ρ[(z.>=2400) .&(z.<2500)] .= 10^0.3
ρ[(z.>=2500) .&(z.<3500)] .= 10^0.7
ρ[(z.>=3500)]              .= 10^1.0
##
transD_GP.CSEM1DInversion.plotmodelfield!(F, z, ρ)
d, σ = transD_GP.addnoise(F, z, ρ, noisefrac=0.05,
               dz=dz, extendfrac=extendfrac, nfixed=nfixed)
savefig("data_model.png", dpi=300)
##
csem = transD_GP.CSEMRadialEr(F, d, σ, z=z, ρ=ρ, nfixed=nfixed)
