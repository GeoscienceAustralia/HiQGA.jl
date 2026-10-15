## do the inversion
m, χ², λ², idx = transD_GP.gradientinv(σstart, σ0, aem, 
                            ;λ²min, λ²max, β², 
                            ntries, nstepsmax, knownvalue,
                            lo, hi, regtype);
