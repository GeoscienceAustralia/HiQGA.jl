using HiQGA.transD_GP
fname_dfn = "L103004_new.dfn"
# electronics file
fname_specs_halt = "electronics_halt.jl"
# column numbers from hdr file
X, Y, Z = "E", "N", "DTM_AHD"
fid = "Fid"
linenum = "Line"
frame_height ="Height"
frame_dx = "dx"
frame_dy = "dy"
frame_dz = "dz"
HM_Z = "HM_Z_"
LM_Z = "LM_Z_"  # first LM gate dropped, also noise and times in electronics_halt
HM_σ = "RUNC_HM_Z_"
LM_σ = "RUNC_LM_Z_"
relerror = false
units = 1e-12
LM_drop,HM_drop = 9,15
##
soundings2 = transD_GP.SkyTEM1DInversion.read_survey_files(fname_dfn,
						             fname_specs_halt = fname_specs_halt,
						             LM_Z             = LM_Z,
									 HM_Z             = HM_Z,
									 LM_σ             = LM_σ,
									 HM_σ             = HM_σ,
									 LM_drop          = LM_drop,
									 HM_drop          = HM_drop,
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