using PyPlot, HiQGA.transD_GP, Random, LinearAlgebra, Statistics

## model
Random.seed!(3)
nlayers, dz = 50, 5.
ρbg        = 1000
zfixed     = [-1e5,   (0:dz:(nlayers-1)*dz )...]
rho        = [1e12,   10 .^(log10(ρbg) .+ randn(nlayers))...]
ianom      = 5
rho[ianom] = 10.
δ          = 1e-4  # fd step size in log10

##  geometry and modeling parameters
rRx = 13.
zRx = -42.0
zTx = -40.0
freqlow = 1e-3
include("electronics_halt.jl")
calcjacobian = true

Fhm = transD_GP.AEM_VMD_HMD.HFieldDHT(
                      lowpassfcs = lowpassfcs,
                      nmax   = length(zfixed),
                      times  = HM_times,
                      ramp   = HM_ramp,
                      zTx    = zTx,
                      rRx    = rRx,
                      rTx    = rTx,
                      zRx    = zRx,
                      freqlow = freqlow,
                      calcjacobian = calcjacobian,
                      )
transD_GP.AEM_VMD_HMD.getfieldTD!(Fhm, zfixed, rho)

## plot HM
f, ax = plt.subplots(1,3, sharey=true, figsize=(8,6))
delz, delt = diff(zfixed[2:end]), diff(log10.(Fhm.times))
tlabel = [(log10.(Fhm.times) .- delt[1]/2)..., log10.(Fhm.times)[end]+delt[end]/2]
zlab = [zfixed[2:end]..., zfixed[end]+delz[end]]
ax[1].step(log10.(rho[2:end]), zfixed[2:end])
im2 = ax[2].pcolormesh(tlabel, zlab, Fhm.dBzdt_J[2:end,:], shading="auto")
cb2 = colorbar(im2, ax=ax[2])
im3 = ax[3].pcolormesh(tlabel, zlab, log10.(abs.(Fhm.dBzdt_J[2:end,:])))
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
J = copy(Fhm.dBzdt_J)

function getJ_TD!(Fin, zin, ρ, Δ)
    # sensitivities computed in log10σ
    # same as in AEM_VMD_HMD code
    Jout = similar(Fin.dBzdt_J)
    modelfd = copy(ρ)
    jflag = Fin.calcjacobian
    Fin.calcjacobian = false
    for i = 2:length(zin)
        modelfd[:] = ρ
        modelfd[i] = 10 .^-(-log10(ρ[i]) .+ Δ/2)
        transD_GP.AEM_VMD_HMD.getfieldTD!(Fin, zin, modelfd)
        HTD_z1 = copy(Fin.dBzdt)
        modelfd[i] = 10 .^-(-log10(ρ[i]) .- Δ/2)
        transD_GP.AEM_VMD_HMD.getfieldTD!(Fin, zin, modelfd)  
        Jout[i,:] = (HTD_z1 - Fin.dBzdt) / Δ
    end
    Fin.calcjacobian = jflag
    Jout
end
J_fd = getJ_TD!(Fhm, zfixed, rho, δ)

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
ax[2].set_xlabel("log10 s")
ax[3].set_xlabel("log10 s")
ax[1].set_ylabel("Depth m")
ax[2].set_title("∂f/∂(log10σ) analytic")
ax[3].set_title("∂f/∂(log10σ) FD")

## more detailed
figure()
nlayer = 40
semilogx(HM_times, J_fd[nlayer,:], "--k")
semilogx(HM_times, J[nlayer,:], "-")

## Benchmark timings
using BenchmarkTools
@info "analytic timings"
@btime transD_GP.AEM_VMD_HMD.getfieldTD!($Fhm, $zfixed, $rho)
@info "Cent Diff timing"
@btime J_fd = getJ_TD!($Fhm, $zfixed, $rho, $δ);