## do the inversion
m, χ², λ², idx = transD_GP.gradientinv(σstart, σ0, tempest, nstepsmax=20, 
                            ;λ²min, λ²max, β², ntries,
                            lo, hi, regtype);
aem = tempest;