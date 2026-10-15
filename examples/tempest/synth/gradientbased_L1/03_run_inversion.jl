## nuisance stuff
if tempest.vectorsum
    nustart  = [zRx-1, x_rx-1]
    nubounds = [zRx-3   zRx+3
                x_rx-3 x_rx+3]
else
    nustart  = [zRx-1, x_rx+1, rx_pitch+0.5]
    nubounds = [zRx-3      zRx+3
                x_rx-3     x_rx+3
                rx_pitch-1.5 rx_pitch+1.5]
end         
## do the inversion
tempest = deepcopy(Torig);
m, nu, χ², χ²nu, λ², idx = transD_GP.gradientinv(σstart, σ0, nustart, tempest; 
                            # Occam stuff
                            λ²min, λ²max, β², ntries,
                            lo, hi, nstepsmax,
                            knownvalue, 
                            L1transform,
                            cd_maxiter=200_000, cd_tol = 1e-6,
                            # optim stuff
                            nubounds, 
                            ntriesnu = 5,
                            boxiters = 2, 
                            );
aem = tempest;