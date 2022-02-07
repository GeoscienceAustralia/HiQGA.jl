## test against Kaufman and Keller 1981 Fig 3.1
using transD_GP.MT1D
## Nomograms
ρ1 = 10.
h1 = 100.
Y = [5000000, 330, 100, 40, 20, 10, 5, 2.5]
ρ2 = [Y.*ρ1;1 ./reverse(Y) .*ρ1]
MT1D.twolayer_ex(h1, ρ1, ρ2)