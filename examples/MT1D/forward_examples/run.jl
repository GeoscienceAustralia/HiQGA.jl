using transD_GP.MT1D, PyPlot
##
ntimesperdecade = 10
T = 10 .^(-3:1/ntimesperdecade:6)
f = 1 ./T
## model 1 conductor above
ρ1 = [10, 100]
z = [0, 5000]
h = diff(z)
Z = MT1D.Z_f(f, ρ1, h)
fig = MT1D.plotcurve(T, Z)
## model 2 conductor below
ρ2 = [100, 10]
Z = MT1D.Z_f(f, ρ2, h)
MT1D.plotcurve(T, Z, fig)
## and to do this with the models plotted
fig = MT1D.plotmodelcurve(T, ρ1, z)
MT1D.plotmodelcurve(T, ρ2, z, fig)