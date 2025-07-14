## read the soundings
cd(@__DIR__)
using PyPlot, HiQGA.transD_GP
fname_dat = "L9000001.XYZ"
zipsaveprefix = splitpath(fname_dat)[end]
fname_specs_halt="electronics_halt.jl"
include(fname_specs_halt)
soundings = transD_GP.TEMPEST1DInversion.read_survey_files(;
	fname_dat,
	times, ramp, Hx_add_noise, Hz_add_noise,
	frame_height = 42,
	frame_dz = 38,
	frame_dx = 37,
	frame_dy = 39,
	Hxs = [1, 15],
	Hzs = [16, 30],
	Hxp = 43,
	Hzp = 44,
	units = 1e-15,
	yaw_rx = 33,
	pitch_rx = 31,
	roll_rx = 32,
	yaw_tx = 36,
	pitch_tx = 34,
	roll_tx = 35,
	X = 40,
	Y = 41,
	Z = 45,
	fid = 46,
	startfrom = 1,
	skipevery = 50,
	linenum = 47,
	makeqcplots=true)

