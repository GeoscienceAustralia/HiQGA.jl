using PyPlot, DelimitedFiles, Random, Statistics,
      HiQGA.transD_GP
## model fixed parts, i.e., air
zfixed   = [-1e5]
ρfixed   = [1e12]
# discretization
zstart = 0.0
extendfrac, dz = 1.03, 1.5
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=65, showplot=false)
z, ρ, nfixed = transD_GP.makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
## geometry and modeling parameters
zTx = -30.0
calcjacobian = true
include("../waveletapprox/electronics_halt.jl")
## make VTEM operator
aem = transD_GP.VTEM1DInversion.dBzdt(;
times, ramp, rTx, zTx, z, ρ, calcjacobian)
## make a model
ρ[(z.>=zstart) .& (z.<50)] .= 20.
ρ[(z.>=50) .&(z.<80)]      .= 1
ρ[(z.>=80) .&(z.<100)]     .= 20
ρ[(z.>=100) .&(z.<200)]    .= 50
ρ[(z.>=200) .&(z.<250)]    .= 80
ρ[(z.>=250)]               .= 150
# add jitter to model in log10 domain
Random.seed!(11)
ρ = 10 .^(0.1*randn(length(ρ)) + log10.(ρ))
# plot it
transD_GP.VTEM1DInversion.plotmodelfield!(aem, log10.(ρ[2:end]))
## make noisy synthetic
transD_GP.VTEM1DInversion.makenoisydata!(aem, log10.(ρ[2:end]); σ_halt)
## Jacobian comparisons
function getJ_TD!(aem_, log10σ, Δ, Jout, modelfd)
    # sensitivities computed in log10σ
    # same as in AEM_VMD_HMD code
    for i = 1:length(log10σ)
        modelfd[:] = log10σ
        modelfd[i] = log10σ[i] .+ Δ/2
        transD_GP.VTEM1DInversion.getfield!(-modelfd, aem_)
        Jout[i,:] = aem_.F.dBzdt[:]
        modelfd[i] = log10σ[i] .- Δ/2
        transD_GP.VTEM1DInversion.getfield!(-modelfd, aem_) 
        Jout[i,:] = (Jout[i,:] - aem_.F.dBzdt) / Δ
    end
end
δ = 1e-2
modelfd = zeros(length(ρ)-1)
aem_ = transD_GP.VTEM1DInversion.dBzdt(;
    times, ramp, rTx, zTx, z, ρ, calcjacobian=false)
J_ = zeros(length(ρ)-1, length(times))
log10σ = -log10.(ρ[2:end])
@btime getJ_TD!($aem_, $log10σ, $δ, $J_, $modelfd)
@btime transD_GP.VTEM1DInversion.getfield!($(-log10σ), $aem)
## plot them
fig, ax = plt.subplots(1,2, sharex=true, sharey=true)
img1 = ax[1].imshow(J_)
colorbar(img1,ax=ax[1])
img2 = ax[2].imshow(aem.J')
colorbar(img2,ax=ax[2])