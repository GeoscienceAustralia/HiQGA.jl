using PyPlot, DelimitedFiles, Random, Statistics, Revise,
      transD_GP
## model fixed parts, i.e., air
ρbg        = 100
zfixed     = [-1e5,   0,   (20:20:200)...]
rho        = [1e12,   ρbg, ρbg*ones(10)...]
ianom      = 5
rho[ianom] = 10.
##  geometry and modeling parameters
rRx = 13.
zRx = -42.0
zTx = -40.0
freqlow = 1e-3
include("electronics_halt.jl")
## LM operator
Flm = transD_GP.AEM_VMD_HMD.HFieldDHT(
                      lowpassfcs = lowpassfcs,
                      times  = LM_times,
                      ramp   = LM_ramp,
                      nmax   = length(zfixed),
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
                      nmax   = length(zfixed),
                      zTx    = zTx,
                      rRx    = rRx,
                      rTx    = rTx,
                      zRx    = zRx,
                      freqlow = freqlow
                      )
## add noise to data
transD_GP.plotmodelfield_skytem!(Flm, Fhm, zfixed, rho)
#remember to specify halt noise if it is in pV
dlow, dhigh, σlow, σhigh = transD_GP.addnoise_skytem(Flm, Fhm,
                zfixed, rho, noisefrac=0.03, halt_LM = LM_noise*1e-12, halt_HM = HM_noise*1e-12,
                dz=diff(zfixed)[1], extendfrac=0.001, nfixed=1)
dlow, dhigh, σlow, σhigh = (dlow, dhigh, σlow, σhigh)./transD_GP.SkyTEM1DInversion.μ₀
aem = transD_GP.dBzdt(Flm, Fhm, dlow, dhigh, σlow, σhigh, z=zfixed, ρ=rho, nfixed=1)
## let's try gradient descent, all model values are in log10 conductivity
σstart, σ0 = map(x->zeros(length(rho)-1), 1:2)
σstart .= -2
σ0 .= -2
λ² = 10 .^range(0, 20, length=10)
## do it
m, χ², idx = transD_GP.gradientinv(σstart, σ0, aem, λ², nstepsmax=10, 
                            regularizeupdate=false,
                            R = transD_GP.makeregR1(aem));
## debug plots
alpha = 0.5
for (i, mi) in enumerate(m)
    figure(figsize=(3,6))
    for (ii, mmi) in enumerate(mi)
        step(mmi, zfixed[2:end], alpha=alpha) 
    end
    step(mi[idx[i]], zfixed[2:end], color="k", linewidth=2)
    step(-log10.(rho[2:end]), zfixed[2:end], color="k", linewidth=2, linestyle="--")
    gca().invert_yaxis()
end