using HiQGA
# If using data below please cite 
# Probabilistic inversion of audio-frequency magnetotelluric 
# data and application to cover thickness estimation for mineral exploration in Australia, JAG
# https://doi.org/10.1016/j.jappgeo.2022.104869
fname = "edifiles/Adv1_AMT.edi"
zfixed   = [-1e5]
ρfixed   = [1e12]
# z grid spec starts, the first z and first ρ will be unused in MT
zstart = 0.0
extendfrac, dz = 1.06, 3.0
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=50, showplot=true)
freqs, d_log10_ρ, d_phase_deg, σ_log10_ρ, σ_phase_deg = transD_GP.MT1DInversion.read_edi(fname, showplot=true, errorfrac=.06)
idx = (10 .<= freqs .<= 1e4) 
freqs, d_log10_ρ, d_phase_deg, σ_log10_ρ, σ_phase_deg = map(x->x[idx],(freqs, d_log10_ρ, d_phase_deg, σ_log10_ρ, σ_phase_deg))
F = transD_GP.MT1DInversion.MT1DZ_nodepthprior(;d_log10_ρ, d_phase_deg, σ_log10_ρ, σ_phase_deg, 
                                    freqs, zboundaries, irxlayer=1, useML=true)
transD_GP.MT1DInversion.plotdata(F)