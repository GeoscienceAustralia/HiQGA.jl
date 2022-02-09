using transD_GP, PyPlot
## model fixed parts, i.e., air, but only to be compatible with AEM ...
zfixed   = [-1e5]
ρfixed   = [1e12]
nmax = 200
# z grid spec starts, the first z and first ρ will be unused in MT
zstart = 0.
extendfrac, dz = 1.028, 10
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=50, showplot=true, atol=1e-3)
z, ρ, nfixed = transD_GP.makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
## fill in detail in ohm-m
ρ[(z.>=zstart) .& (z.<50)] .= 10.
ρ[(z.>=50)  .& (z.<100)] .= 20
ρ[(z.>=100) .& (z.<250)] .= 50
ρ[(z.>=250) .& (z.<300)] .= 5
ρ[(z.>=300) .& (z.<500)] .= 100
ρ[z.>=500] .= 800
##
T = 10 .^range(-3, 0, length=10)
F = transD_GP.MT1DInversion.create_synthetic(ρ[2:end], zboundaries, 1 ./T, rseed=125, showplot=true, noisefrac=0.05, logscaledepth=false)
transD_GP.get_misfit(F.d_log10_ρ, F.d_phase_deg, F.σ_log10_ρ, F.σ_phase_deg, 1 ./T, ρ[2:end], zboundaries)
## now apply a stretch prior
# remember there must be a ρlow and ρhigh at every zall
ρsmooth = log10.(ρ[2:end])
ρsmooth[2:end] = (ρsmooth[1:end-1] + ρsmooth[2:end])/2
ρlow, ρhigh = -0.5*ones(size(zall)), zall*0.001 .+ 2.5
Δ = ρhigh - ρlow
F = transD_GP.MT1DInversion.makestretchop(F, ρlow=ρlow, Δ=Δ)
ax = gcf().axes[1]
transD_GP.MT1DInversion.plotpriorenv(F, ax=ax, lc = "r")
F.stretch = true # if wanting to toggle stretch priors

