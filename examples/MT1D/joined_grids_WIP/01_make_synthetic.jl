using HiQGA.transD_GP, PyPlot
## model fixed parts, i.e., air, but only to be compatible with AEM ...
zfixed   = [-1e5]
ρfixed   = [1e12]
nmax = 200
# z grid spec starts, the first z and first ρ will be unused in MT
zstart = 0.
extendfrac1, dz1, n1 = 1.06, 1.15, 50
# top grid for AEM
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac1, dz=dz1, n=n1, showplot=false, forextension=true)
# deeper grid for MT
extendfrac2, dz2, n2 = 1.5, nothing, 15
zall, znall, zboundaries = transD_GP.zcontinue(;zall=zall, znall=znall, zboundaries=zboundaries, 
                                            extendfrac=extendfrac2, dz=dz2, n=n2, showplot=true)
z, ρ, nfixed = transD_GP.makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
## fill in detail in ohm-m
ρ[(z.>=zstart) .& (z.<50)] .= 20.
ρ[(z.>=50)  .& (z.<100)] .= 10
ρ[(z.>=100) .& (z.<250)] .= 50
ρ[(z.>=250) .& (z.<300)] .= 1
ρ[(z.>=300) .& (z.<500)] .= 80
ρ[(z.>=500) .& (z.<1000)] .= 500
ρ[(z.>=1000) .& (z.<10_000)] .= 700
ρ[z .> 10_000] .= 800
##
T = 10 .^range(-3, 0, length=10)
F = transD_GP.MT1DInversion.create_synthetic(ρ =  ρ[2:end], zboundaries = zboundaries, 
                                            freqs = 1 ./T, rseed=125, showplot=true, noisefrac=0.05, logscaledepth=false)
## now apply a stretch prior
# remember there must be a ρlow and ρhigh at every zall
ρlow, ρhigh = -0.5*ones(size(zall)), zall*0.00005 .+ 2.8
Δ = ρhigh - ρlow
F = transD_GP.MT1DInversion.makestretchop(F, ρlow=ρlow, Δ=Δ)
ax = gcf().axes[1]
transD_GP.MT1DInversion.plotpriorenv(F, ax=ax, lc = "r")
F.stretch = false # if wanting to toggle stretch priors

