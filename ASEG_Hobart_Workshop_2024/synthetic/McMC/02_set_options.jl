# # Prior settings
fileprefix = "SkyTEM_synth_"
# type of Gaussian Process kernel
K = transD_GP.GP.OrstUhn()
# number of GP nuclei
nmin, nmax = 2, 40
# log10 RESISTIVITY bounds to sample between
fbounds = [-1 3.]
# depth locations in number of layers to interpolate resistivity to
xall = copy(permutedims(znall))
# bounds of the depth locations
xbounds = permutedims([extrema(znall)...])
# correlation length in layer number units, tolerance nugget for GP resistivity
λ, δ = [2], 0.1 
# # McMC Proposal specifications as fractions of prior bounds
sdev_pos = 0.05*vec(diff([extrema(znall)...]))
sdev_prop = [0.07*diff(fbounds, dims=2)...]
sdev_dc = [0.01*diff(fbounds, dims=2)...]
##
# # If restarting
restart = false
if !restart # delete earlier runs if starting afresh
    history_mode = "w"
    deletefiles = ["misfits_", "points_", "models_"].*fileprefix.*"s.bin"
    [isfile(f) && rm(f) for f in deletefiles]
else
    history_mode = "a"
end
##
# # Initialize a stationary GP using these options
opt = transD_GP.OptionsStat(;
            dispstatstoscreen = true, # show/don't stats in Jupyter
            nmin, nmax, xbounds, fbounds, fdataname = fileprefix,
            xall, λ, δ, sdev_prop, sdev_pos, sdev_dc, history_mode,
            quasimultid = false, save_freq = 50, K);