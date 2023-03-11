## make options for the lengthscale GP
# Prior stuff
log10bounds = [-1.1 -0.9] # length scale bounds in log10 of distance
δlog10λ = 0.02 # nugget of lengthscales in distnce units
nminlog10λ, nmaxlog10λ = 2, 20 # no. of GP nuclei for lengthscales
Klog10λ = transD_GP.GP.OrstUhn() # kernel for lengthscales
λlog10λ = [0.02abs(diff([extrema(x)...])[1])] # lengthscale for lengthscales, is stationary
demeanλ = false # must be centred about lengthscale of lengthscales -- cannot demean or sample this
sampledcλ = false # must be centred about lengthscale of lengthscales -- cannot demean or sample this
xall = permutedims(collect(x)) # discretization
xbounds = permutedims([extrema(x)...]) # discretization bounds
# McMC proposal stuff
sdev_poslog10λ = [0.05abs(diff([extrema(x)...])[1])]
sdev_proplog10λ = 0.05*diff(log10bounds, dims=2)[:]
## Initialize a lengthscale model using these options
Random.seed!(12)
optlog10λ = transD_GP.OptionsStat(nmin = nminlog10λ,
                        nmax = nmaxlog10λ,
                        xbounds = xbounds,
                        fbounds = log10bounds,
                        xall = xall,
                        λ = λlog10λ,
                        δ = δlog10λ,
                        demean = demeanλ,
                        sampledc = sampledcλ,
                        sdev_prop = sdev_proplog10λ,
                        sdev_pos = sdev_poslog10λ,
                        needλ²fromlog = true, # needed for variable lengthscales
                        updatenonstat = true, # needed for variable lengthscales
                        quasimultid = false,
                        K = Klog10λ,
                        save_freq = 20,
                        timesλ = 4,
                        peskycholesky = true
                        )
## make options for the nonstationary actual properties GP
# prior stuff
nmin, nmax = 2, 30 # no. of nuclei
fbounds = permutedims([extrema(ynoisy[.!isnan.(ynoisy)])...]) # property bounds
ymin, ymax = extrema(y) 
fbounds[1] > ymin && (fbounds[1] = ymin)
fbounds[2] < ymax && (fbounds[2] = ymax)
# proposal stuff
sdev_prop = 0.05*diff(fbounds, dims=2)[:]
sdev_pos = [0.02abs(diff([extrema(x)...])[1])]
sdev_dc = 0.008*diff(fbounds, dims=2)[:]
K = transD_GP.GP.Mat32() # Kernel for nonstationary GP
δ = 0.05 # property nugget
opt = transD_GP.OptionsNonstat(optlog10λ;
                        nmin = nmin,
                        nmax = nmax,
                        fbounds = fbounds,
                        δ = δ,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        K = K,
                        sdev_dc
                        )