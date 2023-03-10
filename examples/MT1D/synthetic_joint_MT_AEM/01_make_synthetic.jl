using HiQGA.transD_GP, PyPlot, Random, Statistics
## model fixed parts, i.e., air, but only to be compatible with AEM ...
zfixed   = [-1e5]
ρfixed   = [1e12]
# z grid spec starts, the first z and first ρ will be unused in MT
zstart = 0.
extendfrac, dz = 1.03, 1.5
nlayers = 100
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=nlayers, showplot=true)
z, ρ, nfixed = transD_GP.makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
## fill in detail in ohm-m
ρ[(z.>=zstart) .& (z.<50)] .= 20.
ρ[(z.>=50) .&(z.<80)]      .= 1
ρ[(z.>=80) .&(z.<100)]     .= 20
ρ[(z.>=100) .&(z.<200)]    .= 50
ρ[(z.>=200) .&(z.<250)]    .= 80
ρ[(z.>=250) .& (z.<300)]   .= 2
ρ[(z.>=300) .& (z.<400)]   .= 800
ρ[(z.>=400) .& (z.<500)]   .= 500
ρ[(z.>=500) .& (z.<600)]   .= 80
ρ[(z.>=600) .& (z.<700)]   .= 400
ρ[z.>=700]                 .= 800
# add jitter
Random.seed!(42)
ρ = 10 .^(0.1*randn(length(ρ)) + log10.(ρ))
## AMT operator
freqs = 10 .^range(log10(11.2), log10(8800), length=37)
# plot with resistivity axis reversed 
FMT = transD_GP.MT1DInversion.create_synthetic(;ρ =  ρ[2:end],zboundaries, revax=true,
                                            freqs, showplot=true, noisefrac=0.05, logscaledepth=false)
## VTEM operator
zTx = -30.0
include("../../VTEM/synth/waveletapprox/electronics_halt.jl")
aem = transD_GP.VTEM1DInversion.dBzdt(;times, ramp, rTx, zTx, z, ρ)
transD_GP.VTEM1DInversion.makenoisydata!(aem, log10.(ρ[2:end]); σ_halt)
## Joint Operator
F = transD_GP.GenericJointInversion.JointOperator([aem, FMT])