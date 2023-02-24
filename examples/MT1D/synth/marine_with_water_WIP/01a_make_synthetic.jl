using HiQGA.transD_GP, PyPlot
## model fixed parts, i.e., air, but only to be compatible with AEM ...
zfixed   = [-1e5]
ρfixed   = [1e12]
nmax = 200
# z grid spec starts, the first z and first ρ will be unused in MT
zstart = 0.
extendfrac, dz = 1.028, 10
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=50, showplot=true)
z, ρ, nfixed = transD_GP.makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
## fill in detail in ohm-m
ρwater = 0.22
ρ[(z.>=zstart) .& (z.<35)] .= ρwater
ρ[(z.>=35) .& (z.<100)] .= 4.
ρ[(z.>=100) .& (z.<250)] .= 14
ρ[(z.>=250) .&(z.<600)] .= 0.8
ρ[z.>=600] .= 100
##
# ATTENTION: to get the first non water layer if wanting to use with CSEM
irxlayer = findlast(isequal.(ρ[2:end], ρwater))+1
T = 10 .^range(-2, 0, length=10)
F = transD_GP.MT1DInversion.create_synthetic(ρ[2:end], zboundaries, 1 ./T, irxlayer=irxlayer,
        rseed=125, showplot=true, noisefrac=0.05, logscaledepth=false)
transD_GP.get_misfit(F.d_log10_ρ, F.d_phase_deg, F.σ_log10_ρ, F.σ_phase_deg, 1 ./T, ρ[2:end], zboundaries, irxlayer)