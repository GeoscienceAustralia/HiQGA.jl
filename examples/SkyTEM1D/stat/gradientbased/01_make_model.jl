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
regtype = :R1
## do it
m, χ², λ², idx = transD_GP.gradientinv(σstart, σ0, aem, nstepsmax=15, 
                            regularizeupdate=false, 
                            dobo=true, 
                            λ²min=0,
                            λ²max=7, 
                            ntries=8,
                            knownvalue=0.7, 
                            frac=5, 
                            breakonknown = true,
                            firstvalue=:last, 
                            κ = transD_GP.GP.Mat52(),
                            regtype = regtype);
## debug plots: all in each
alpha = 0.1
idxlast=length(m)
for (i, mi) in enumerate(m)
    if isempty(λ²[i]) 
        idxlast = i-1
        break
    end    
    f = figure(figsize=(7,6))
    s1 = subplot(131)
    for ii in 1:length(χ²[i])
        s1.step(-mi[ii], aem.z[2:end], alpha=alpha) 
    end
    s1.step(-mi[idx[i]], aem.z[2:end], color="r", linewidth=2)
    s1.step(log10.(rho[2:end]), aem.z[2:end], color="k", linewidth=2, linestyle="--")
    s1.set_ylabel("depth m")
    s1.set_xlabel("Log₁₀ ρ")
    s1.invert_yaxis()
    s1.invert_xaxis()
    s2 = subplot(232)
    hps = hcat(λ²[i]...)'
    sortedλidx = sortperm(hps[:,1])
    sortedαidx = sortperm(hps[:,2])
    s2.loglog(hps[sortedλidx,1], χ²[i][sortedλidx], ".-", markersize=5)
    s2.plot(hps[idx[i],1], χ²[i][idx[i]], "o" )
    s2.set_ylabel("χ²")
    s3 = subplot(236)
    s3.loglog(χ²[i][sortedαidx], hps[sortedαidx,2], ".-", markersize=5)
    s3.plot(χ²[i][idx[i]], hps[idx[i],2], "o" )
    s3.set_xlabel("χ²")
    s4 = subplot(235, sharex=s2, sharey=s3)
    s4.scatter(hps[:,1], hps[:,2], c=log10.(χ²[i]), cmap="magma")
    s4.plot(hps[idx[i],1], hps[idx[i],2], "x", markersize=10)
    s4.set_xlabel("λ²"); s4.set_ylabel("step length α")
    plt.tight_layout()
end
## debug plots best in each
figure(figsize=(3,6))
alpha = 0.1
for (i, mi) in enumerate(m)
    step(-mi[idx[i]], aem.z[2:end], alpha=alpha)
    i == idxlast && break
end
step(-m[idxlast][idx[idxlast]], aem.z[2:end], color="r", linewidth=2)
step(log10.(rho[2:end]), aem.z[2:end], color="k", linewidth=2, linestyle="--")
gca().invert_yaxis()
gca().invert_xaxis()
ylabel("depth m")
xlabel("Log₁₀ ρ")
plt.tight_layout()
ax = gca()
## Compare with posterior model covariance
using LinearAlgebra
F = aem
R = transD_GP.makereg(regtype, F)
JtW, Wr = F.J'*F.W, F.W*F.res
H = JtW*(JtW)' + λ²[idxlast][idx[idxlast]]*R'R
Cpost = inv(Hermitian(H))
figure()
zplot = [aem.z[2:end]; aem.z[end] + 5]
pcolormesh(zplot, zplot, Cpost)
xlabel("Depth m")
ylabel("Depth m")
colorbar()
sd = sqrt.(diag(Cpost))
ax.fill_betweenx(zplot[1:end-1]+diff(zplot)/2, -m[idxlast][idx[idxlast]] -sd, -m[idxlast][idx[idxlast]] +sd, alpha=0.2)
ax.set_xlim(3, -1)