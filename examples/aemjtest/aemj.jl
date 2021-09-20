using PyPlot, transD_GP, Random
## frequencies
nFreqsPerDecade     = 5
freqLowLimit        = 1e-3
freqHighLimit       = 1e5
freqs = 10 .^(log10(freqLowLimit):1/nFreqsPerDecade:log10(freqHighLimit))
## model
zfixed   = [-1e5,   0,   (20:20:200)...]
rho      = [1e12,   100, 100*ones(10)...]
rho[8] = 1.
rho[5] = 1.
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