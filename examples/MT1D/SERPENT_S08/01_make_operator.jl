using transD_GP, PyPlot, DelimitedFiles, DataInterpolations
## model fixed parts, i.e., air, but only to be compatible with AEM ...
zfixed   = [-1e5]
ρfixed   = [1e12]
nmax = 200
# z grid spec starts, the first z and first ρ will be unused in MT
zstart = 0.0
extendfrac, dz = 1.169, 10
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=50, showplot=true, atol=1e-3)
z, ρ, nfixed = transD_GP.makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
## SERPENT S08
D = readdlm("MTS08.txt")
T = 1 ./D[:,1]
d_log10_ρ = D[:,2]
d_phase_deg = -D[:,3]
σ_log10_ρ = D[:,4]
σ_phase_deg = D[:,5]
## load zrho limits
usedepthprior = false # toggle this and see differences in posterior
zrholim = readdlm("zlog10rholims.txt")
interplow, interphigh = map(x->ConstantInterpolation(x, zrholim[:,1]), (zrholim[:,2], zrholim[:,3]))
ρlow, ρhigh = interplow.(zall), interphigh.(zall)
Δ = ρhigh - ρlow
F = transD_GP.MT1DInversion.MT1DZ_nodepthprior(d_log10_ρ, d_phase_deg,σ_log10_ρ, σ_phase_deg, 1 ./T, zboundaries, 1)
F = transD_GP.MT1DInversion.makestretchop(F, ρlow=ρlow, Δ=Δ) # apply the stretch prior and get new operator
transD_GP.MT1DInversion.plotdata(F)
F.stretch = usedepthprior 
if usedepthprior
    transD_GP.MT1DInversion.plotpriorenv(F, lc = "r", plotlinear=false)
end    