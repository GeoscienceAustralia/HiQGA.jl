## stationary GP options
# McMC prior stuff
xall = transD_GP.get_image_prediction_points(img) # image discretization
xbounds = [img.x[1] img.x[end];img.y[1] img.y[end]] # image discretization bounds
nmin, nmax= 2, 100 # no. of GP nuclei to represent pixel values with
K = transD_GP.GP.Mat32() # Kernel function for GP
ftrain = img.d[img.select] # training pixel values
fbounds = permutedims([extrema(ftrain)...]) # bounds of pixel values
δ = 0.3 # pixel nugget
λ = [141, 141] # lengthscales for pixel GP
# McMC proposal stuff
sdev_prop = 0.05*diff(fbounds, dims=2)[:]
sdev_pos = 0.05*abs.([diff([extrema(xall[1,:])...])[1], diff([extrema(xall[2,:])...])[1]])
sdev_dc = 0.01*diff(fbounds, dims=2)[:]
## Initialize a stationary GP using these options
opt = transD_GP.OptionsStat(;nmin = nmin,
                        nmax = nmax,
                        xbounds,
                        fbounds = fbounds,
                        xall = xall,
                        λ = λ,
                        δ = δ,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        sdev_dc,
                        quasimultid = false,
                        K = K,
                        )
