using HiQGA
cd(@__DIR__)
# common to both
zfixed   = [-1e5]
ρfixed   = [1e12]
zstart = 0.0
extendfrac, dz = 1.03, 1.5
nlayers = 100
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=nlayers, showplot=true)
## read MT data
fname = "WIT2.edi"
freqs, d_log10_ρ, d_phase_deg, σ_log10_ρ, σ_phase_deg = transD_GP.MT1DInversion.read_edi(fname, showplot=true, errorfrac=.06)
idx = (10 .<= freqs .<= 1e4) 
freqs, d_log10_ρ, d_phase_deg, σ_log10_ρ, σ_phase_deg = map(x->x[idx],(freqs, d_log10_ρ, d_phase_deg, σ_log10_ρ, σ_phase_deg))
FMT = transD_GP.MT1DInversion.MT1DZ_nodepthprior(;d_log10_ρ, d_phase_deg, σ_log10_ρ, σ_phase_deg, 
                                    freqs, zboundaries, irxlayer=1, useML=false)
transD_GP.MT1DInversion.plotdata(FMT)
## info to read VTEM data
fname_dat = "ST_Line_1131.dat"
zipsaveprefix = splitpath(fname_dat)[end]
# electronics file
fname_specs_halt = "electronics_halt.jl"
# column numbers from hdr file
X, Y, Z = 28, 29, 31
fid = 5
linenum = 4
frame_height = 30
d = [177,221]
##
using HiQGA.transD_GP.NearestNeighbors
include(fname_specs_halt)
soundings = transD_GP.VTEM1DInversion.read_survey_files(; X, Y, Z, 
									fid, linenum, frame_height,
									d, fname_dat, lowpassfcs, times, 
									ramp, σ_halt, rTx,
									startfrom        = 1,
									skipevery        = 1,
									dotillsounding   = nothing,
									makeqcplots      = true)
X_MT, Y_MT = 275876.48, 6842372.28								
sno = transD_GP.CommonToAll.getclosestidx(X_MT, Y_MT, soundings)
soundings = soundings[sno:sno]
