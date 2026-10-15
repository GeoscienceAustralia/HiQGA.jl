## for gradient descent, all model values are in log10 conductivity
# create empty log10 conductivity arrays for start and background
σstart      = -2*ones(length(zall)) 
σ0          = -2*ones(length(zall)) # prior model
lo          = -3      # log10 min conductivity  
hi          = 1       # log10 max conductivity
λ²min       = 1       # min log10 of Tikhonov parameter
λ²max       = 3.      # max log10 of Tikhonov parameter
ntries      = 10      # how many max tries between λ²min, λ²max
β²          = 0.001   # what fraction of roughness penalty to use to enforce a penalty in deviations from background
knownvalue  = 0.7     # cutoff to start next Gauss-Newton step
nstepsmax   = 20      # Gauss-Newton steps
L1transform = transD_GP.R1Transform(aem) # L1 regularisation