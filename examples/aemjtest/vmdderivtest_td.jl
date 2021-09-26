using PyPlot, transD_GP, Random, LinearAlgebra

## model
ρbg        = 1000
zfixed     = [-1e5,   0,   (20:20:200)...]
rho        = [1e12,   ρbg, ρbg*ones(10)...]
ianom      = 5
rho[ianom] = 10.
δ          = 1e-4  # fd step size in log10

## Gauss Newton settings
use_fd          = false # in the one step for Gauss Newton
startfromcloser = false  # whether to start midway between anomaly and bg or plain bg in G-N step
λ²              = 5e-17 # model update damping

##  geometry and frequencies for VMD in frequency domain
nFreqsPerDecade     = 5
freqLowLimit        = 1e-3
freqHighLimit       = 1e5
freqs = 10 .^(log10(freqLowLimit):1/nFreqsPerDecade:log10(freqHighLimit))
rRx = 100.
zRx = -37.
zTx = -35.
nkᵣeval = 50
modelprimary = true
doconvramp = false
provideddt = true
Fvmd = transD_GP.AEM_VMD_HMD.HFieldDHT(
                      calcjacobian = true,
                      nmax = length(zfixed),
                      zTx    = zTx,
                      rRx    = rRx,
                      freqs  = freqs,
                      zRx    = zRx,
                      nkᵣeval = nkᵣeval,
                      modelprimary = modelprimary,
                      doconvramp = doconvramp,
                      provideddt = provideddt)
transD_GP.AEM_VMD_HMD.getfieldTD!(Fvmd, zfixed, rho)

## plot
times = Fvmd.interptimes
f, ax = plt.subplots(1,3, sharey=true, figsize=(8,6))
delz, delt = diff(zfixed[2:end]), diff(log10.(freqs))
tlabel = [(log10.(times) .- delt[1]/2)..., log10.(times)[end]+delt[end]/2]
zlab = [zfixed[2:end]..., zfixed[end]+delz[end]]
ax[1].step(log10.(rho[2:end]), zfixed[2:end])
im2 = ax[2].pcolormesh(tlabel, zlab, Fvmd.HTD_z_J_interp[2:end,:])
cb2 = colorbar(im2, ax=ax[2])
im3 = ax[3].pcolormesh(tlabel, zlab, log10.(abs.(Fvmd.HTD_z_J_interp[2:end,:])))
cb2 = colorbar(im3, ax=ax[3])
ax[1].set_xlabel("log10 Ω-m")
ax[1].invert_xaxis()
ax[2].invert_yaxis()
ax[2].set_xlabel("log10 s")
ax[3].set_xlabel("log10 s")
ax[1].set_ylabel("Depth m")
ax[2].set_title("∂f/∂(log10σ)")
ax[3].set_title("∂f/∂(log10σ)")

## approximate partials with central differences and compare
J = copy(Fvmd.HTD_z_J_interp)

function getJ_TD_noconv!(Fin, zin, ρ, Δ)
    # sensitivities computed in log10σ
    # same as in AEM_VMD_HMD code
    Jout = similar(Fin.HTD_z_J_interp)
    for i = 2:length(zin)
        modelfd = copy(ρ)
        modelfd[i] = 10 .^-(-log10(ρ[i]) .+ Δ/2)
        transD_GP.AEM_VMD_HMD.getfieldTD!(Fin, zin, modelfd)
        HTD_z1 = copy(Fin.HTD_z_interp)
        modelfd[i] = 10 .^-(-log10(ρ[i]) .- Δ/2)
        transD_GP.AEM_VMD_HMD.getfieldTD!(Fin, zin, modelfd)  
        Jout[i,:] = (HTD_z1 - Fin.HTD_z_interp) / Δ
        end
    Jout
end
J_fd = getJ_TD_noconv!(Fvmd, zfixed, rho, δ)

## comparison plots
f, ax = plt.subplots(1,3, sharey=true, figsize=(10,6))
ax[1].step(log10.(rho[2:end]), zfixed[2:end])
im2 = ax[2].pcolormesh(tlabel, zlab, J[2:end,:])
cb2 = colorbar(im2, ax=ax[2])
im3 = ax[3].pcolormesh(tlabel, zlab, J_fd[2:end,:])
cb2 = colorbar(im3, ax=ax[3])
ax[1].set_xlabel("log10 Ω-m")
ax[1].invert_xaxis()
ax[2].invert_yaxis()
ax[2].set_xlabel("log10 Hz")
ax[3].set_xlabel("log10 s")
ax[1].set_ylabel("Depth m")
ax[2].set_title("∂f/∂(log10σ) analytic")
ax[3].set_title("∂f/∂(log10σ) FD")
