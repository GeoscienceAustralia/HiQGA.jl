using transD_GP, PyPlot
## model fixed parts, i.e., air, but only to be compatible with AEM ...
zfixed   = [-1e5]
ρfixed   = [1e12]
nmax = 200
# z grid spec starts, the first z and first ρ will be unused in MT
zstart = 0.0
extendfrac, dz = 1.22, 15
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=50, showplot=true, atol=1e-3)
z, ρ, nfixed = transD_GP.makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
## COPROD from Occam paper
T = [28.5 38.5 52.0 70.5 95.5 129.0 174.6 236.2 319.6 432.5 585.1 791.7 1071.1 1449.2 1960.7][:]
d_log10_ρ = [2.315  2.254  2.229  2.188  2.18  2.162  2.151  2.208  2.194  2.299  2.338  2.42  2.405  2.308  2.397][:]
d_phase_deg = -[57.19  58.19  61.39  59.09  59.89  51.19  46.89  42.79  36.89  32.0  44.0  32.0  37.59  45.29  50.09][:]
σ_log10_ρ = [0.0721  0.0425  0.0244  0.021  0.0164  0.0173  0.0287  0.0328  0.0193  0.027  0.0591  0.0506  0.0825  0.1233  0.0927][:]
σ_phase_deg = [22.95  22.95  4.96  4.46  5.96  22.95  22.95  2.46  1.65  22.95  6.37  2.46  22.95  4.15  22.95][:]
F = transD_GP.MT1DInversion.MT1DZ_nodepthprior(d_log10_ρ, d_phase_deg,σ_log10_ρ, σ_phase_deg, 1 ./T, zboundaries, 1)
transD_GP.MT1DInversion.plotdata(F)