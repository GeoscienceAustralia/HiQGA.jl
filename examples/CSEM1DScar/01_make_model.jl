srcdir = dirname(dirname(pwd()))*"/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
using PyPlot, DelimitedFiles, Random, MAT, Revise,
      CSEM1DEr, GeophysOperator
##
fdataname = "s26.mat"
file = matopen(fdataname)
freqs   = read(file, "f")[:] # Hz
rRx   = read(file, "r")[:]  # Ranges to receivers   (m)
zRx = read(file, "zRx")[1]     # Depth of Receivers, if this is an array it calls normal FHT without lagged convolution
zTx = read(file, "zTx")
RxAzim = read(file, "RxAz")
TxDip = read(file, "TxDip")
zfixed   = read(file, "z")[:]
ρfixed = read(file, "rho")[:]
σ = read(file, "sd")
noisyd = read(file, "data")
abovesflr = read(file, "abovesflr")

# Note that the receiver depth needs to be in same model layer as transmitter.
zstart = zRx + abovesflr
extendfrac, dz = 1.008, 12.
zall, znall, zboundaries = GeophysOperator.setupz(zstart, extendfrac, dz=dz, n=126)
z, ρ, nfixed = GeophysOperator.makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
useML = true
##
F = CSEM1DEr.RadialErLagged(zTx    = [zTx],
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
GeophysOperator.plotmodelfield!(F, z, ρ)
##
GeophysOperator.plotmodelfield!(F::CSEM1DEr.RadialEr, z, ρ, noisyd, σ,
                                dz=dz, extendfrac=extendfrac, nfixed=nfixed)
csem = GeophysOperator.CSEMRadialEr(F, noisyd, σ, z=z, ρ=ρ, nfixed=nfixed, useML=useML)
