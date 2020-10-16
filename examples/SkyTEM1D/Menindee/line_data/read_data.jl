srcdir = dirname(dirname(dirname(dirname(pwd()))))*"/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
## column numbers from hdr file
fname_dat = "em_north_line.dat"
frame_height = 9
frame_dx = 21
frame_dy = 22
frame_dz = 23
LM_Z = [26,42]
HM_Z = [43,67]
# electronics file
fname_specs_halt = "electronics_halt.jl"
##
using Revise, SkyTEM1DInversion
SkyTEM1DInversion.read_survey_files(fname_dat         = fname_dat,
						             fname_specs_halt = fname_specs_halt,
						             LM_Z             = LM_Z,
									 HM_Z             = HM_Z,
									 frame_height     = frame_height,
									 frame_dz         = frame_dz,
									 frame_dy         = frame_dy,
									 frame_dx         = frame_dx)
##
ss = SkyTEM1DInversion.read_survey_files(fname_dat         = fname_dat,
						             fname_specs_halt = fname_specs_halt,
						             LM_Z             = LM_Z,
									 HM_Z             = HM_Z,
									 frame_height     = frame_height,
									 frame_dz         = frame_dz,
									 frame_dy         = frame_dy,
									 frame_dx         = frame_dx,
									 dotillsounding=2,
									 makesounding=true)
