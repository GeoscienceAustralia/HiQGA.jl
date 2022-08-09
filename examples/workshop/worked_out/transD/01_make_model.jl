using PyPlot, DelimitedFiles, Random, Statistics, Printf
using HiQGA.transD_GP
## model description
# fixed parts, i.e., air
zfixed   = [-1e5]
ρfixed   = [1e12]
# part to be inverted
zstart = 0.0 # surface height, i.e., ground
dz = 1.5 # first cell thickness
extendfrac = 1.03 # extend each cell by this much with depth
# make a 65 layer model
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, showplot=true, n=65)
z, ρ, nfixed = transD_GP.makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
##  geometry and modeling parameters
rRx = 13. # radial distance to Rx 
zRx = -42.0 # distance above ground to Rx 
zTx = -40.0 # distance above ground to Tx
calcjacobian = false # turn OFF for transD otherwise very slow. Must be on for gradient-based
include("../../exercises/synth/gradientbased/electronics_halt.jl") # waveforms, gates, etc.
## LM operator
Flm = transD_GP.AEM_VMD_HMD.HFieldDHT(
                      lowpassfcs = lowpassfcs,
                      times  = LM_times,
                      ramp   = LM_ramp,
                      zTx    = zTx,
                      rRx    = rRx,
                      rTx    = rTx,
                      zRx    = zRx,
                      calcjacobian = calcjacobian
                      )
## HM operator
Fhm = transD_GP.AEM_VMD_HMD.HFieldDHT(
                      lowpassfcs = lowpassfcs,
                      times  = HM_times,
                      ramp   = HM_ramp,
                      zTx    = zTx,
                      rRx    = rRx,
                      rTx    = rTx,
                      zRx    = zRx,
                      calcjacobian = calcjacobian
                      )
## define model and plot the forward response
ρ[(z.>=zstart) .& (z.<50)] .= 20.
ρ[(z.>=50) .&(z.<80)] .= 1
ρ[(z.>=80) .&(z.<100)] .= 20
ρ[(z.>=100) .&(z.<200)] .= 50
ρ[(z.>=200) .&(z.<250)] .= 80
ρ[(z.>=250)]            .= 150
# add jitter to model in log10 domain
Random.seed!(11)
ρ = 10 .^(0.1*randn(length(ρ)) + log10.(ρ))
# plot the forward calculation
transD_GP.plotmodelfield_skytem!(Flm, Fhm, z, ρ)
## add noise to data
# remember to specify halt noise if it is in pV
dlow, dhigh, σlow, σhigh = transD_GP.addnoise_skytem(Flm, Fhm,
    z, ρ, noisefrac=0.03, halt_LM = LM_noise*1e-12, halt_HM = HM_noise*1e-12,
    dz=dz, extendfrac=extendfrac, nfixed=nfixed)
# EM data for inversion are always in units of H = B/μ
dlow, dhigh, σlow, σhigh = (dlow, dhigh, σlow, σhigh)./transD_GP.SkyTEM1DInversion.μ₀
## create inversion operator with physics, data, noise and model discretization 
aem = transD_GP.dBzdt(Flm, Fhm, dlow, dhigh, σlow, σhigh, z=z, ρ=ρ, nfixed=nfixed);
