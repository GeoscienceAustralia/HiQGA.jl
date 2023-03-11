## lengthscale GP options
# Prior stuff
xall = transD_GP.get_image_prediction_points(img) # image discretization
nminlog10λ, nmaxlog10λ = 2, 100 # no. of lengthscale muclei
Klog10λ = transD_GP.GP.Mat32() # GP kernel type for lengthscales
log10bounds = [log10(10) log10(2000); log10(10) log10(2000)] # bounds for vaiable lengthscales in log 10 units
δlog10λ = 0.3 # nugget of lengthscales in log 10 units
λlog10λ = 0.1*abs.([diff([extrema(xall[1,:])...])[1], diff([extrema(xall[2,:])...])[1]]) # stationary lengthscale for lengthscales
demean = false # essential for length scale change approximations
sampledc = false # can't have jumping dc value for lengthscale approximation
# McMC proposal stuff
sdev_proplog10λ = 0.05*diff(log10bounds, dims=2)[:]
sdev_poslog10λ = 0.05*abs.([diff([extrema(xall[1,:])...])[1], diff([extrema(xall[2,:])...])[1]])
## Initialize a lengthscale model using these options
optlog10λ = transD_GP.OptionsStat(;nmin = nminlog10λ,
                        nmax = nmaxlog10λ,
                        xbounds = [img.x[1] img.x[end];img.y[1] img.y[end]],
                        fbounds = log10bounds,
                        xall = xall,
                        λ = λlog10λ,
                        δ = δlog10λ,
                        needλ²fromlog = true, # needed for variable lengthscales
                        updatenonstat = true, # needed for variable lengthscales
                        demean, sampledc,
                        sdev_prop = sdev_proplog10λ,
                        sdev_pos = sdev_poslog10λ,
                        quasimultid = false,
                        K = Klog10λ,
                        timesλ = 3.6
                        )
## make options for the nonstationary actual properties GP
nmin, nmax = 2, 100 # no. of nuclei for property (i.e., pixel) values
ftrain = img.d[img.select] # data to regress
fbounds = permutedims([extrema(ftrain)...]) # bounds of properties (i.e. pixels)
δ = 0.25 # pixel value nugget
K = transD_GP.GP.Mat32() # Pixel values GP kernel function
# McMC proposal stuff
sdev_prop = 0.05*diff(fbounds, dims=2)[:] 
sdev_pos = 0.05*abs.([diff([extrema(xall[1,:])...])[1], diff([extrema(xall[2,:])...])[1]])
sdev_dc = 0.01*diff(fbounds, dims=2)[:]
## Initialize model for the nonstationary properties GP
Random.seed!(13)
opt = transD_GP.OptionsNonstat(optlog10λ;
                        nmin = nmin,
                        nmax = nmax,
                        fbounds = fbounds,
                        δ = δ,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        sdev_dc,
                        K = K
                        )
