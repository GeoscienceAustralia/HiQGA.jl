using HiQGA, DelimitedFiles
zstart = 0.0
extendfrac, dz = 1.06, 1.25
nlayers = 52
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=nlayers, showplot=true)
zcwell, rhowell = map(i->readdlm("testing_data/BHMAR23-1.con", skipstart=6)[:,i], 1:2)
rhowell = 3 .-log10.(rhowell) # log10 ohm m
##
Mwell = transD_GP.CommonToAll.block1Dvalues(
        [rhowell], zcwell, [zboundaries[1:end-1] zboundaries[2:end]], :mean
)
