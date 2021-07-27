## info to read data
# datafile
fname_dat = "L103004.XYZ"
# electronics file
fname_specs_halt = "electronics_halt.jl"
# column numbers from hdr file
X, Y, Z = 2, 70, 1
fid = 3
linenum = 43
frame_height = 4
frame_dx = 71
frame_dy = 72
frame_dz = 73
HM_Z = [20, 42]
LM_Z = [52+1, 69] # first LM gate dropped, also noise and times in electronics_halt
relerror = false
units = 1e-12
##
using transD_GP
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
									 units            = units,
									 relerror         = relerror,
									 linenum          = linenum,
                                     startfrom        = 1,
									 skipevery        = 34,
									 dotillsounding   = nothing,
									 makesounding     = true)
