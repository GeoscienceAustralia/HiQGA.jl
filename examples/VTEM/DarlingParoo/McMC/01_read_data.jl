using HiQGA.transD_GP
## info to read data
# datafile
fname_dat = "../gradientbased/200.asc"
zipsaveprefix = splitpath(fname_dat)[end]
# electronics file
fname_specs_halt = "../gradientbased/electronics_halt.jl"
# column numbers from hdr file
X, Y, Z = 28, 29, 634
fid = 5
linenum = 4
frame_height = 30
d = [177,221]
##
soundings = transD_GP.VTEM1DInversion.read_survey_files(; X, Y, Z, 
									fid, linenum, frame_height,
									d, fname_dat, fname_specs_halt,
									tx_rx_dz_pass_through = 0.01, # GA-AEM Z up, +ve is rx over tx
									startfrom        = 1,
									skipevery        = 1,
									dotillsounding   = 3,
									makeqcplots      = true)
