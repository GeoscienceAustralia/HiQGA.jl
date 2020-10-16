srcdir = dirname(dirname(dirname(dirname(pwd()))))*"/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
## info to read data
# datafile
fname_dat = "em_north_line.dat"
# electronics file
fname_specs_halt = "electronics_halt.jl"
# column numbers from hdr file
X, Y = 12, 13
fid = 1
linenum = 2
frame_height = 9
frame_dx = 21
frame_dy = 22
frame_dz = 23
LM_Z = [26,42]
HM_Z = [43,67]
##
using Revise, SkyTEM1DInversion
SkyTEM1DInversion.read_survey_files(fname_dat         = fname_dat,
						             fname_specs_halt = fname_specs_halt,
						             LM_Z             = LM_Z,
									 HM_Z             = HM_Z,
									 frame_height     = frame_height,
									 frame_dz         = frame_dz,
									 frame_dy         = frame_dy,
									 frame_dx         = frame_dx,
									 X                = X,
									 Y                = Y,
									 fid              = fid,
									 linenum          = linenum)
##
ss = SkyTEM1DInversion.read_survey_files(fname_dat         = fname_dat,
						             fname_specs_halt = fname_specs_halt,
						             LM_Z             = LM_Z,
									 HM_Z             = HM_Z,
									 frame_height     = frame_height,
									 frame_dz         = frame_dz,
									 frame_dy         = frame_dy,
									 frame_dx         = frame_dx,
									 X                = X,
									 Y                = Y,
									 fid              = fid,
									 linenum          = linenum,
									 skipevery        = 2,
									 dotillsounding   = 4,
									 makesounding     = true)
