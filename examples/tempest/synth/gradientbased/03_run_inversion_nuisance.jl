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
                            nstepsmax=30, 
                            # Occam stuff
                            λ²min, λ²max, β², ntries,
                            lo, hi, regtype,
                            # optim stuff
                            nubounds, 
                            ntriesnu = 5,
                            boxiters = 2, 
                            usebox = true,
                            reducenuto = 0.2,
                            debuglevel = 2,
                            breaknuonknown = false);
aem = tempest;
## for plotting forwards
mn = reduce(hcat, [transD_GP.TEMPEST1DInversion.setnuforinvtype(tempest, n) for n in nu])'
mfinal = [-mm[i] for (mm, i) in zip(m, idx)]
mtest = copy(mfinal[end])
nutest = copy(mn[end,:])
if tempest.vectorsum
    @info transD_GP.TEMPEST1DInversion.get_misfit(-mtest, [zRx, x_rx], tempest, nubounds)
else
    @info transD_GP.TEMPEST1DInversion.get_misfit(-mtest, [zRx, x_rx, rx_pitch], tempest, nubounds)
end        
@info transD_GP.TEMPEST1DInversion.get_misfit(-mtest, nu[end], tempest, nubounds)
transD_GP.plotmodelfield!(tempest, mtest, nutest)