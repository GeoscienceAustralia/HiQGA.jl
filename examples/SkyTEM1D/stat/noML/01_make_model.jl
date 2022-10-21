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
nmax = length(z)
## LM operator
Flm = transD_GP.AEM_VMD_HMD.HFieldDHT(
                      lowpassfcs = lowpassfcs,
                      times  = LM_times,
                      ramp   = LM_ramp,
                      nmax   = nmax,
                      zTx    = zTx,
                      rRx    = rRx,
                      rTx    = rTx,
                      zRx    = zRx,
                      freqlow = freqlow
                      )
## HM operator
Fhm = transD_GP.AEM_VMD_HMD.HFieldDHT(
                      lowpassfcs = lowpassfcs,
                      times  = HM_times,
                      ramp   = HM_ramp,
                      nmax   = nmax,
                      zTx    = zTx,
                      rRx    = rRx,
                      rTx    = rTx,
                      zRx    = zRx,
                      freqlow = freqlow
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
## add noise to data
transD_GP.plotmodelfield_skytem!(Flm, Fhm, z, ρ)
#remember to specify halt noise if it is in pV
dlow, dhigh, σlow, σhigh = transD_GP.addnoise_skytem(Flm, Fhm,
                z, ρ, noisefrac=0.03, halt_LM = LM_noise*1e-12, halt_HM = HM_noise*1e-12,
                dz=dz, extendfrac=extendfrac, nfixed=nfixed)
## create operator
# but with a coarser grid
extendfrac, dz = 1.06, 1.15
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=50, showplot=true)
zgrid, ρgrid, nfixed = transD_GP.makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
dlow, dhigh, σlow, σhigh = (dlow, dhigh, σlow, σhigh)./transD_GP.SkyTEM1DInversion.μ₀
aem = transD_GP.dBzdt(Flm, Fhm, dlow, dhigh, σlow, σhigh, z=zgrid, ρ=ρgrid, nfixed=nfixed)
