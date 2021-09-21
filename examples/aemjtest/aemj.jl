using PyPlot, transD_GP, Random
## frequencies
nFreqsPerDecade     = 5
freqLowLimit        = 1e-3
freqHighLimit       = 1e5
freqs = 10 .^(log10(freqLowLimit):1/nFreqsPerDecade:log10(freqHighLimit))
## model
zfixed   = [-1e5,   0,   (20:20:200)...]
rho      = [1e12,   1000, 1000*ones(10)...]
rho[5] = 1000.
rho[10] = 1000.
##  geometry
rRx = 100.
zRx = -37.
zTx = -35.
nkᵣeval = 50
## VMD
modelprimary = true
Fvmd = transD_GP.AEM_VMD_HMD.HFieldDHT(
                      calcjacobian = true,
                      nmax = length(zfixed),
                      zTx    = zTx,
                      rRx    = rRx,
                      freqs  = freqs,
                      zRx    = zRx,
                      nkᵣeval = nkᵣeval,
                      modelprimary = modelprimary)
transD_GP.AEM_VMD_HMD.getfieldFD!(Fvmd, zfixed, rho)

## plot
f, ax = plt.subplots(1,3, sharey=true, figsize=(8,6))
delz, delf = diff(zfixed[2:end]), diff(log10.(freqs))
flabel = [(log10.(freqs) .- delf[1]/2)..., log10.(freqs)[end]+delf[end]/2]
zlab = [zfixed[2:end]..., zfixed[end]+delz[end]]
ax[1].step(log10.(rho[2:end]), zfixed[2:end])
im2 = ax[2].pcolormesh(flabel, zlab, abs.(Fvmd.HFD_z_J[2:end,:]), shading="auto")
cb2 = colorbar(im2, ax=ax[2])
im3 = ax[3].pcolormesh(flabel, zlab, log10.(abs.(Fvmd.HFD_z_J[2:end,:])))
cb2 = colorbar(im3, ax=ax[3])
ax[1].set_xlabel("log10 Ω-m")
ax[1].invert_xaxis()
ax[2].invert_yaxis()
ax[2].set_xlabel("log10 Hz")
ax[3].set_xlabel("log10 Hz")
ax[1].set_ylabel("Depth m")
ax[2].set_title("∂f/∂(log10σ)")
ax[3].set_title("∂f/∂(log10σ)")

## approximate partials with central differences:
fd_pertubation_cond = 1e-6
J_fd = similar(Fvmd.HFD_z_J)
J = copy(Fvmd.HFD_z_J)
HFD_z = copy(Fvmd.HFD_z)

@time for i = 2:length(zfixed)
    modelfd = copy(rho)
    modelfd[i] = 10 .^-(-log10(rho[i]) .+ fd_pertubation_cond/2)
    transD_GP.AEM_VMD_HMD.getfieldFD!(Fvmd, zfixed, modelfd)
    HFD_z1 = copy(Fvmd.HFD_z)
    modelfd[i] = 10 .^-(-log10(rho[i]) .- fd_pertubation_cond/2)
    transD_GP.AEM_VMD_HMD.getfieldFD!(Fvmd, zfixed, modelfd)  
    J_fd[i,:] = (HFD_z1 - Fvmd.HFD_z) / fd_pertubation_cond
end

## comparison plots
f, ax = plt.subplots(1,3, sharey=true, figsize=(10,6))
ax[1].step(log10.(rho[2:end]), zfixed[2:end])
im2 = ax[2].pcolormesh(flabel, zlab, abs.(J[2:end,:]))
cb2 = colorbar(im2, ax=ax[2])
im3 = ax[3].pcolormesh(flabel, zlab, abs.(J_fd[2:end,:]))
cb2 = colorbar(im3, ax=ax[3])
ax[1].set_xlabel("log10 Ω-m")
ax[1].invert_xaxis()
ax[2].invert_yaxis()
ax[2].set_xlabel("log10 Hz")
ax[3].set_xlabel("log10 Hz")
ax[1].set_ylabel("Depth m")
ax[2].set_title("∂f/∂(log10σ) analytic")
ax[3].set_title("∂f/∂(log10σ) FD")