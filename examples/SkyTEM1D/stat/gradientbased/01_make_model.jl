using PyPlot, DelimitedFiles, Random, Statistics, Revise,
      transD_GP
## model fixed parts, i.e., air
ρbg        = 100
zfixed     = [-1e5,   0,   (20:10:200)...]
rho        = [1e12,   ρbg*ones(length(zfixed)-1)...]
ianom      = 5
rho[ianom] = 10
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
σstart .= -2.5
σ0 .= -2
λ² = 10 .^range(0, 8, length=20)
## do it
m, χ², idx = transD_GP.gradientinv(σstart, σ0, aem, λ², nstepsmax=10, 
                            regularizeupdate=false,
                            regtype = :R0);
## debug plots: all in each
alpha = 0.5
for (i, mi) in enumerate(m)
    figure(figsize=(7,6))
    for (ii, mmi) in enumerate(mi)
        subplot(121)
        step(-mmi, zfixed[2:end], alpha=alpha) 
    end
    step(-mi[idx[i]], zfixed[2:end], color="r", linewidth=2)
    step(log10.(rho[2:end]), zfixed[2:end], color="k", linewidth=2, linestyle="--")
    ylabel("depth m")
    xlabel("Log₁₀ ρ")
    gca().invert_yaxis()
    gca().invert_xaxis()
    subplot(122)
    loglog(λ², χ²[i])
    plot(λ²[idx[i]], χ²[i][idx[i]], "." )
    plt.tight_layout()
end
## debug plots best in each
figure(figsize=(3,6))
alpha = 0.25
for (i, mi) in enumerate(m)
    step(-mi[idx[i]], zfixed[2:end], alpha=alpha)
end
step(-m[end][idx[end]], zfixed[2:end], color="r", linewidth=2)
step(log10.(rho[2:end]), zfixed[2:end], color="k", linewidth=2, linestyle="--")
gca().invert_yaxis()
gca().invert_xaxis()
ylabel("depth m")
xlabel("Log₁₀ ρ")
plt.tight_layout()