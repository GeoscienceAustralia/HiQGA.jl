## BO nuisance stuff
nustart  = [zRx-2, x_rx+2]
nuλ²frac = [8, 8]
nubounds = [zRx-2.5   zRx+2.5
            x_rx-2.5 x_rx+2.5]
ndivsnu  = [30, 30]
ntriesnu = 8 
## do the inversion
tempest = deepcopy(Torig);
m, nu, χ², χ²nu, λ², idx, idxnu = transD_GP.gradientinv(σstart, σ0, nustart, tempest, nstepsmax=20, 
                            ;nubounds, nuλ²frac, ndivsnu, ntriesnu,
                            λ²min, λ²max, β², ntries,
                            lo, hi, regtype);
aem = tempest;