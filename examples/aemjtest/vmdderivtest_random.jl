using PyPlot, transD_GP, Random, LinearAlgebra, Statistics

## model
Random.seed!(3)
nlayers, dz = 50, 5.
ρbg        = 1000
zfixed     = [-1e5,   (0:dz:(nlayers-1)*dz )...]
rho        = [1e12,   10 .^(log10(ρbg) .+ randn(nlayers))...]
ianom      = 5
rho[ianom] = 10.
δ          = 1e-4  # fd step size in log10

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

## approximate partials with central differences and compare
J = copy(Fvmd.HFD_z_J)

function getJ!(Fin, zin, ρ, Δ)
    # sensitivities computed in log10σ
    # same as in AEM_VMD_HMD code
    Jout = similar(Fin.HFD_z_J)
    modelfd = copy(ρ)
    jflag = Fin.calcjacobian
    Fin.calcjacobian = false
    for i = 2:length(zin)
        modelfd[:] = ρ
        modelfd[i] = 10 .^-(-log10(ρ[i]) .+ Δ/2)
        transD_GP.AEM_VMD_HMD.getfieldFD!(Fin, zin, modelfd)
        HFD_z1 = copy(Fin.HFD_z)
        modelfd[i] = 10 .^-(-log10(ρ[i]) .- Δ/2)
        transD_GP.AEM_VMD_HMD.getfieldFD!(Fin, zin, modelfd)  
        Jout[i,:] = (HFD_z1 - Fin.HFD_z) / Δ
    end
    Fin.calcjacobian = jflag
    Jout
end
J_fd = getJ!(Fvmd, zfixed, rho, δ)

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

## more detailed
figure()
nlayer = 3
subplot(121)
semilogx(freqs, real(J_fd[nlayer,:]), "--k")
semilogx(freqs, real(J[nlayer,:]), "-")
subplot(122)
semilogx(freqs, imag(J_fd[nlayer,:]), "--k")
semilogx(freqs, imag(J[nlayer,:]), "-")