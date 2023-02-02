## for gradient descent, all model values are in log10 conductivity
# create empty log10 conductivity arrays for start and background
σstart, σ0 = map(x->zeros(length(tempest.ρ)-nfixed), 1:2)
σstart .= -2. # constant value assigned to multiple elements of array (start)
σ0 .= -2 # constant value assigned to multiple elements of array (background)
regtype = :R1 # :R1 is first difference regularization, :R2 is two succesive first differences
lo, hi = -3, 1 # min, max log10 conductivity
λ²min, λ²max = -0.5, 7 # min, max of regularization parameter (Tikhonov parameter)
ntries = 8 # how many max tries between λ²min, λ²max
β² = 0. # what fraction of roughness penalty to use to enforce a penalty in deviations from background