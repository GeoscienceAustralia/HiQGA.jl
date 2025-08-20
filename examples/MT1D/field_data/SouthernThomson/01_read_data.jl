using HiQGA
cd(@__DIR__)
## read MT data
fname = "WIT2.edi"
freqs, d_log10_ρ, d_phase_deg, σ_log10_ρ, σ_phase_deg = transD_GP.MT1DInversion.read_edi(fname, showplot=true, errorfrac=.06)
idx = (10 .<= freqs .<= 1e4) 
freqs, d_log10_ρ, d_phase_deg, σ_log10_ρ, σ_phase_deg = map(x->x[idx],(freqs, d_log10_ρ, d_phase_deg, σ_log10_ρ, σ_phase_deg))
## make z grid
zfixed   = [-1e5]
ρfixed   = [1e12] # first value not used
zstart = 0.0
extendfrac, dz = 1.03, 1.5
nlayers = 100
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=nlayers, showplot=true)
FMT = transD_GP.MT1DInversion.MT1DZ_nodepthprior(;d_log10_ρ, d_phase_deg, σ_log10_ρ, σ_phase_deg, 
                                    freqs, zboundaries, irxlayer=1, useML=false)
transD_GP.MT1DInversion.plotdata(FMT)