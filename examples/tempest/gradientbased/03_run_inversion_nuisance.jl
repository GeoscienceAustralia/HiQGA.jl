## BO nuisance stuff
if tempest.vectorsum
    nustart  = [zRx-1, x_rx+1]
    nuλ²frac = [8, 8]
    nubounds = [zRx-2.5   zRx+2.5
                x_rx-2.5 x_rx+2.5]
    ndivsnu  = [30, 30]
    ntriesnu = 8
else
    nustart  = [zRx-2, x_rx+2, rx_pitch+0.5]
    nuλ²frac = [8, 8, 8]
    nubounds = [zRx-2.5      zRx+2.5
                x_rx-2.5     x_rx+2.5
                rx_pitch-0.5 rx_pitch+0.5]
    ndivsnu  = [30, 30, 14]
    ntriesnu = 10
end         
## do the inversion
tempest = deepcopy(Torig);
m, nu, χ², χ²nu, λ², idx = transD_GP.gradientinv(σstart, σ0, nustart, tempest, nstepsmax=30, 
                            ;nubounds, nuλ²frac, ndivsnu,
                            λ²min, λ²max, β², ntries,
                            # optim stuff
                            ntriesnu, 
                            breaknuonknown=true,
                            reducenuto=0.2,
                            lo, hi, regtype);
aem = tempest;
## for plotting forwards
mn = reduce(hcat, [transD_GP.TEMPEST1DInversion.setnuforinvtype(tempest, n) for n in nu])'
mfinal = [-mm[i] for (mm, i) in zip(m, idx)]
mtest = copy(mfinal[end])
nutest = copy(mn[end,:])
@info transD_GP.TEMPEST1DInversion.get_misfit(-mtest, [zRx, x_rx], tempest, nubounds)
@info transD_GP.TEMPEST1DInversion.get_misfit(-mtest, nu[end], tempest, nubounds)
transD_GP.TEMPEST1DInversion.plotmodelfield!(tempest, mtest, nutest)