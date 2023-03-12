## lengthscale GP options
# McMC prior stuff
xall = transD_GP.get_image_prediction_points(img) # image discretization
nminλ, nmaxλ = 2, 300 # no. of lengthscale muclei
boundsλ = log10.([20. 1000;20  1000]) # bounds for vaiable lengthscales in log 10 units
Kλ = transD_GP.GP.OrstUhn() # GP kernel type for lengthscales
δλ = diff(img.x)[1] # nugget of lengthscales in log 10 units
λλ = 0.1*abs.([diff([extrema(xall[1,:])...])[1], diff([extrema(xall[2,:])...])[1]]) # stationary lengthscale for lengthscales
demean = false # essential for length scale change approximations
sampledc = false # can't have jumping dc value for lengthscale approximation
# McMC proposal stuff
sdev_propλ = 0.05*diff(boundsλ, dims=2)[:]
sdev_posλ = 0.05*abs.([diff([extrema(xall[1,:])...])[1], diff([extrema(xall[2,:])...])[1]])
## Initialize a lengthscale model using these options
optλ = transD_GP.OptionsStat(;nmin = nminλ,
                        nmax = nmaxλ,
                        xbounds = [minimum(xall[1,:]) maximum(xall[1,:]);
                                   minimum(xall[2,:]) maximum(xall[2,:])],
                        fbounds = boundsλ,
                        xall = xall,
                        λ = λλ,
                        δ = δλ,
                        needλ²fromlog = true, # needed for variable lengthscales
                        updatenonstat = true, # needed for variable lengthscales
                        demean, sampledc,
                        sdev_prop = sdev_propλ,
                        sdev_pos = sdev_posλ,
                        quasimultid = false,
                        K = Kλ,
                        timesλ = 3.8,
                        history_mode = "w" # a for append, "w" when starting fresh
                        )
## make options for the nonstationary actual properties GP
nmin, nmax = 2, 300 # no. of nuclei for property (i.e., pixel) values
ftrain = img.d[img.select] # data to regress
fbounds = permutedims([extrema(ftrain)...]) # bounds of properties (i.e. pixels)
δ = 0.1 # pixel value nugget
K = transD_GP.GP.OrstUhn() # Pixel values GP kernel function
# McMC proposal stuff
sdev_prop = 0.05*diff(fbounds, dims=2)[:]
sdev_pos = 0.05*abs.([diff([extrema(xall[1,:])...])[1], diff([extrema(xall[2,:])...])[1]])
sdev_dc = 0.01*diff(fbounds, dims=2)[:]
## Initialize model for the nonstationary properties GP
opt = transD_GP.OptionsNonstat(optλ;
                        nmin = nmin,
                        nmax = nmax,
                        fbounds = fbounds,
                        δ = δ,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        sdev_dc,
                        K = K
                        )
