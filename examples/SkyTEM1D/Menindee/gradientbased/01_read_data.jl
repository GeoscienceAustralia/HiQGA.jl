## info to read data
# datafile
fname_dat = "../McMC/coincidewithtempest.dat"
zipsaveprefix = splitpath(fname_dat)[end]
# electronics file
fname_specs_halt = "../McMC/electronics_halt.jl"
# column numbers from hdr file
X, Y, Z = 12, 13, 15
fid = 1
linenum = 2
frame_height = 9
frame_dx = 21
frame_dy = 22
frame_dz = 23
LM_Z = [26,42]
HM_Z = [43,67]
relerror = false
units = 1e-12
##
using HiQGA.transD_GP
soundings = transD_GP.SkyTEM1DInversion.read_survey_files(fname_dat = fname_dat,
						             fname_specs_halt = fname_specs_halt,
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
									 linenum          = linenum,
                                     startfrom        = 1,
									 skipevery        = 50,
									 dotillsounding   = nothing,
									 relerror         = relerror,
									 units            = units,     
									 makesounding     = true)