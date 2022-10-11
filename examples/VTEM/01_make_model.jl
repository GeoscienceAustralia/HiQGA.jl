using PyPlot, DelimitedFiles, Random, Statistics, Revise,
      HiQGA.transD_GP
## model fixed parts, i.e., air
zfixed   = [-1e5]
ρfixed   = [1e12]
nmax = 200
# Note that the receiver and transmitter need to be in layer 1
zstart = 0.0
extendfrac, dz = 1.03, 1.5
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=65, showplot=true)
z, ρ, nfixed = transD_GP.makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
##  geometry and modeling parameters
rRx = 0.0
zRx = -30.01
zTx = -30.0
freqlow = 1e-3
calcjacobian = true
include("electronics_halt.jl")
## LM operator
F = transD_GP.AEM_VMD_HMD.HFieldDHT(;
                      times,
                      ramp,
                      nmax,
                      zTx,
                      rRx,
                      rTx,
                      zRx,
                      freqlow,
                      calcjacobian
                      )
## fill in detail in ohm-m
ρ[(z.>=zstart) .& (z.<50)] .= 20.
ρ[(z.>=50) .&(z.<80)] .= 1
ρ[(z.>=80) .&(z.<100)] .= 20
ρ[(z.>=100) .&(z.<200)] .= 50
ρ[(z.>=200) .&(z.<250)] .= 80
ρ[(z.>=250)]            .= 150
# add jitter to model in log10 domain
Random.seed!(11)
ρ = 10 .^(0.1*randn(length(ρ)) + log10.(ρ))
# get field
@time transD_GP.AEM_VMD_HMD.getfieldTD!(F, z, ρ)