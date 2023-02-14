## do the inversion
m, χ², λ², idx = transD_GP.gradientinv(σstart, σ0, aem, nstepsmax=20, 
                            ;λ²min, λ²max, β², ntries,
                            lo, hi, regtype);
## if using bo
m, χ², λ², idx = transD_GP.gradientinv(σstart, σ0, aem, nstepsmax=20, 
                            ;β², ntries,
                            lo, hi, regtype,
                            # BO stuff
                            λ²min, λ²max, 
                            αmin = -3, 
                            αmax = 0, 
                            αfrac = 4,
                            ntestdivsα = 32,
                            knownvalue = 0.7,
                            firstvalue = :last,
                            κ = transD_GP.GP.Mat32(),
                            breakonknown = true,
                            dobo = true);

