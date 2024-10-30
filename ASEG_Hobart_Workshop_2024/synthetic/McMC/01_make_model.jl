using PyPlot, HiQGA, Random
cd(@__DIR__)
##
# # Geometry
# model fixed parts, i.e., air
Random.seed!(23)
zfixed   = [-1e5]
ρfixed   = [1e12]
##
# # Model discretization
# Note that the receiver and transmitter need to be in layer 1
zstart = 0.0
extendfrac, dz = 1.06, 1.5
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=50, showplot=true)
z, ρ, nfixed = transD_GP.makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed);
rRx = 13.
##
# # AEM modeling parameters
# Rx-Tx geometry 
zRx = -42.0
zTx = -40.0
include("/scratch/ns59/HiQGA.jl/ASEG_Hobart_Workshop_2024/UDF_data/electronics_halt.jl")
calcjacobian = false # switch off for McMC!
# make SkyTEM operator
aem = transD_GP.SkyTEM1DInversion.dBzdt(;
    timeslow = LM_times, ramplow = LM_ramp, zRxlow=zRx, zTxlow = zTx,
    timeshigh = HM_times, ramphigh = HM_ramp, zRxhigh=zRx, zTxhigh = zTx,
    rRx, rTx, z, ρ, lowpassfcs, calcjacobian);
##
# # Create resistivity model in ohm-m
ρ[(z.>=zstart) .& (z.<50)] .= 20.
ρ[(z.>=50)      .&(z.<80)] .= 1
ρ[(z.>=80)     .&(z.<120)] .= 20
ρ[(z.>=120)    .&(z.<150)] .= 1
ρ[(z.>=150)    .&(z.<250)] .= 50
ρ[(z.>=250)]               .= 150
# add jitter to model in log10 domain
Random.seed!(10)
ρ = 10 .^(0.1*randn(length(ρ)) + log10.(ρ))
##
# # Plot model and data
# plot noise free data due to model
transD_GP.SkyTEM1DInversion.plotmodelfield!(aem, log10.(ρ[2:end]))
# add noise to data
transD_GP.SkyTEM1DInversion.makenoisydata!(aem, log10.(ρ[2:end]); 
    σ_halt_low=LM_noise, σ_halt_high=HM_noise,
    units = 1e-12)
