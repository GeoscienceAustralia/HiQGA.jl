using PyPlot, DelimitedFiles, Random, Statistics, Revise,
      HiQGA.transD_GP
## model fixed parts, i.e., air
zfixed   = [-1e5]
ρfixed   = [1e12]
# discretization
zstart = 0.0
extendfrac, dz = 1.03, 1.5
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=65, showplot=false)
z, ρ, nfixed = transD_GP.makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
## geometry and modeling parameters
zTx = -30.0
calcjacobian = true
include("../electronics_halt.jl")
## make LM operator
aem = transD_GP.VTEM1DInversion.dBzdt(;
times, ramp, rTx, zTx, z, ρ)
## make a model
ρ[(z.>=zstart) .& (z.<50)] .= 20.
ρ[(z.>=50) .&(z.<80)]      .= 1
ρ[(z.>=80) .&(z.<100)]     .= 20
ρ[(z.>=100) .&(z.<200)]    .= 50
ρ[(z.>=200) .&(z.<250)]    .= 80
ρ[(z.>=250)]               .= 150
# add jitter to model in log10 domain
Random.seed!(11)
ρ = 10 .^(0.1*randn(length(ρ)) + log10.(ρ))
## make noisy synthetic
transD_GP.VTEM1DInversion.makenoisydata!(aem, log10.(ρ[2:end]); σ_halt)
   