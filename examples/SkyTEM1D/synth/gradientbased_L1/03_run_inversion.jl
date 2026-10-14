## do the inversion
m, χ², λ², idx = transD_GP.gradientinv(σstart, σ0, aem; 
                            nstepsmax, λ²min, λ²max, β², ntries, 
                            knownvalue, 
                            L1transform,
                            cd_maxiter=200_000, cd_tol=1e-6,
                            lo, hi);