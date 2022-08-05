using PyPlot, DelimitedFiles, Random, Statistics, Revise, Printf,
HiQGA.transD_GP
## model fixed parts, i.e., air
Random.seed!(23)
zfixed   = [-1e5]
ρfixed   = [1e12]
nmax = 200
# Note that the receiver and transmitter need to be in layer 1
zstart = 0.0
extendfrac, dz = 1.03, 1.5
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=65)
z, ρ, nfixed = transD_GP.makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
##  geometry and modeling parameters
rRx = 13.
zRx = -42.0
zTx = -40.0
freqlow = 1e-3
include("electronics_halt.jl")
calcjacobian = true
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
                      freqlow = freqlow,
                      calcjacobian = calcjacobian
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
                      freqlow = freqlow,
                      calcjacobian = calcjacobian
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
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=50)
zgrid, ρgrid, nfixed = transD_GP.makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
dlow, dhigh, σlow, σhigh = (dlow, dhigh, σlow, σhigh)./transD_GP.SkyTEM1DInversion.μ₀
aem = transD_GP.dBzdt(Flm, Fhm, dlow, dhigh, σlow, σhigh, z=zgrid, ρ=ρgrid, nfixed=nfixed)
## let's try gradient descent, all model values are in log10 conductivity
σstart, σ0 = map(x->zeros(length(aem.ρ)-1), 1:2)
σstart .= -2.
σ0 .= -2
regtype = :R1
lo, hi = -3, 1
λ²min, λ²max=-0.5, 7
β² = 0.0
ntries = 6
## do it
m, χ², λ², idx = transD_GP.gradientinv(σstart, σ0, aem, nstepsmax=20, 
                            ;λ²min, λ²max, β², ntries,
                            lo, hi, regtype)
## debug plots: all in each
alpha = 0.3
idxlast=length(m)β²
ndata = length(aem.res)
prefix = "regtype_$(regtype)_β²_$(β²)_iteration_"
for (i, mi) in enumerate(m)
    if isempty(λ²[i]) 
        global idxlast = i-1
        break
    end    
    f, s = plt.subplots(1, 2, gridspec_kw=Dict("width_ratios" => [1,1.5]),
        figsize=(8,6))
    for ii in 1:length(χ²[i])
        s[1].step(-mi[ii], aem.z[2:end], alpha=alpha) 
    end
    s[1].step(-mi[idx[i]], aem.z[2:end], color="r", linewidth=2)
    s[1].step(log10.(ρ[2:end]), z[2:end], color="k", linewidth=2, linestyle="--")
    s[1].set_ylabel("depth m")
    s[1].set_xlabel("Log₁₀ ρ")
    s[1].invert_yaxis()
    s[1].invert_xaxis()
    hps = hcat(λ²[i]...)'
    sortedλidx = sortperm(hps[:,1])
    sortedαidx = sortperm(hps[:,2])
    s[2].loglog(hps[sortedλidx,1], χ²[i][sortedλidx], ".-", markersize=5)
    s[2].plot(hps[idx[i],1], χ²[i][idx[i]], "o" )
    s[2].text(.95hps[idx[i],1], 1.05χ²[i][idx[i]], "$(@sprintf("%.2f", χ²[i][idx[i]]))" )
    s[2].plot(10 .^[λ²min, λ²max], [ndata, ndata], "--k")
    s[2].set_ylabel("χ²")
    s[2].set_xlabel("λ²")
    plt.suptitle("Iteration $i, α=$(hps[idx[i],2])")
    transD_GP.nicenup(gcf(), fsize=12)
    savefig(prefix*"$(@sprintf("%02i",i)).png", dpi=300)
end
## debug plots best in each
figure(figsize=(3,6))
alpha = 0.1
for (i, mi) in enumerate(m)
    step(-mi[idx[i]], aem.z[2:end], alpha=alpha)
    i == idxlast && break
end
step(-m[idxlast][idx[idxlast]], aem.z[2:end], color="r", linewidth=2)
step(log10.(ρ[2:end]), z[2:end], color="k", linewidth=2, linestyle="--")
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
H = JtW*(JtW)' + λ²[idxlast][idx[idxlast]][1]*R'R
Cpost = inv(Hermitian(H))
figure()
zplot = [zboundaries; zboundaries[end] + 5]
pcolormesh(zplot, zplot, Cpost)
xlabel("Depth m")
ylabel("Depth m")
colorbar()
sd = sqrt.(diag(Cpost))
ax.fill_betweenx(zall, -m[idxlast][idx[idxlast]] -sd, -m[idxlast][idx[idxlast]] +sd, alpha=0.2)
ax.set_xlim(3, -1)