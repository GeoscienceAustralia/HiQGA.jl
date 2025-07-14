## info to read data
cd(@__DIR__)
# datafile
fname_dat = "/scratch/ns59/HiQGA.jl/ASEG_Hobart_Workshop_2024/UDF_data/twolines.dat"
# electronics file
fname_specs_halt = "/scratch/ns59/HiQGA.jl/ASEG_Hobart_Workshop_2024/UDF_data/electronics_halt.jl"
# column numbers from hdr file
X, Y, Z = 1, 2, 3
fid = 73
linenum = 72
frame_height = 4
frame_dx = 5
frame_dy = 6
frame_dz = 7
LM_Z = [16, 33]
HM_Z = [49, 71] 
relerror = false
units = 1e-12
##
using HiQGA.transD_GP
include(fname_specs_halt)
soundings = transD_GP.SkyTEM1DInversion.read_survey_files(;
									 fname_dat = fname_dat,
						             LM_Z             = LM_Z,
									 HM_Z             = HM_Z,
									 frame_height     = frame_height,
									 frame_dz         = frame_dz,
									 frame_dy         = frame_dy,
									 frame_dx         = frame_dx,
									 X                = X,
									 Y                = Y,
                                     Z                = Z,
									 fid              = fid,
									 units            = units,
									 relerror         = relerror,
									 linenum          = linenum,
                                     startfrom        = 1,
									 skipevery        = 1,
									 dotillsounding   = nothing,
									 makeqcplots      = true,
									 lowpassfcs,
    								 LM_times, HM_times, LM_ramp, HM_ramp, 
									 LM_noise, HM_noise, rTx)
