using PyPlot, DelimitedFiles, Random, Statistics, Revise,
      HiQGA.transD_GP
## model fixed parts, i.e., air
Random.seed!(23)
zfixed   = [-1e5]
ρfixed   = [1e12]
# Note that the receiver and transmitter need to be in layer 1
zstart = 0.0
extendfrac, dz = 1.03, 1.5
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=65, showplot=true)
z, ρ, nfixed = transD_GP.makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
##  geometry and modeling parameters
rRx = 13.
zRx = -42.0
zTx = -40.0
freqlow = 1e-3
include("electronics_halt.jl")
calcjacobian = true
## make SkyTEM operator
aem = transD_GP.SkyTEM1DInversion.dBzdt(;
    timeslow = LM_times, ramplow = LM_ramp, zRxlow=zRx, zTxlow = zTx,
    timeshigh = HM_times, ramphigh = HM_ramp, zRxhigh=zRx, zTxhigh = zTx,
    rRx, rTx, z, ρ, lowpassfcs, calcjacobian)
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
# plot it
transD_GP.SkyTEM1DInversion.plotmodelfield!(aem, log10.(ρ[2:end]))
## add noise to data
transD_GP.SkyTEM1DInversion.makenoisydata!(aem, log10.(ρ[2:end]); 
    σ_halt_low=LM_noise, σ_halt_high=HM_noise)
