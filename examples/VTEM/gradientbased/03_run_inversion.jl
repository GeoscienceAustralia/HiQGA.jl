## do the inversion
m, χ², λ², idx = transD_GP.gradientinv(σstart, σ0, aem, nstepsmax=20, 
                            ;λ²min, λ²max, β², ntries,
                            lo, hi, regtype);