using PyPlot, Test, Random, Revise, Statistics, LinearAlgebra
any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
import GP, TransD_GP
## make options for the multichannel lengthscale GP
nminlog10λ, nmaxlog10λ = 2, 20
pnorm = 2.
Klog10λ = GP.Mat32()
λx = 1
x = 0:(0.01λx):λx
λlog10λ = [0.05maximum(x)]
demean = false
sdev_poslog10λ = [0.01maximum(x)]
log10bounds = [-1 -0.69]
δlog10λ = 0.1
sdev_proplog10λ = [0.1]
xall = permutedims(collect(x))
xbounds = permutedims([extrema(x)...])
## Initialize a model using these options
Random.seed!(12)
optlog10λ = TransD_GP.Options(nmin = nminlog10λ,
                        nmax = nmaxlog10λ,
                        xbounds = xbounds,
                        fbounds = log10bounds,
                        xall = xall,
                        λ = λlog10λ,
                        δ = δlog10λ,
                        demean = demean,
                        sdev_prop = sdev_proplog10λ,
                        sdev_pos = sdev_poslog10λ,
                        pnorm = pnorm,
                        quasimultid = false,
                        K = Klog10λ,
                        timesλ = 2.
                        )
@time  log10λ = TransD_GP.init(optlog10λ)
## make options for the nonstationary GP
nmin, nmax = 2, 20
fbounds = [-2. 2]
δ = 0.1
sdev_prop = [0.1]
sdev_pos = [0.05]
demean_ns = true
K = GP.Mat32()
## Initialize model for the nonstationary GP
Random.seed!(13)
opt = TransD_GP.Options(nmin = nmin,
                        nmax = nmax,
                        xbounds = optlog10λ.xbounds,
                        fbounds = fbounds,
                        xall = optlog10λ.xall,
                        λ = λlog10λ,
                        δ = δ,
                        demean = demean_ns,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        pnorm = pnorm,
                        quasimultid = false,
                        K = K
                        )
mns = TransD_GP.init(opt, log10λ)
## update λ and mns.fstar through birth in λ
doall = false
lold = copy(log10λ.fstar)
TransD_GP.birth!(log10λ, optlog10λ, mns, opt, doall=doall)
lsc = TransD_GP.testupdate(optlog10λ, log10λ)
ftest = TransD_GP.testupdate(opt, log10λ, mns)

f, ax = plt.subplots(2, 1, sharex=true)
ax[1].plot(xall', sqrt.(log10λ.fstar)')
ax[1].plot(xall', sqrt.(lold)')
ax[1].plot(xall', 10 .^lsc, "--k")
ax[1].plot(log10λ.xtrain[:,1:log10λ.n]', 10 .^(log10λ.ftrain[1:log10λ.n]), "+m")
ax[1].plot(log10λ.xtrain[:,log10λ.n]', 10^(log10λ.ftrain[log10λ.n]), ".r", markersize=10)
ax[1].set_title("λ")

ax[2].plot(xall', mns.fstar)
ax[2].plot(xall', ftest)
ax[2].plot(mns.xtrain[:,1:mns.n]', mns.ftrain[1:mns.n], "xr")
ax[2].set_title("function")

ax2 = ax[1].twinx()
ax2.plot(xall', log10.(abs.(lold' -10 .^2lsc)))
ax[1].grid()
ax[2].grid()
## update mns.fstar through birth in mns
TransD_GP.birth!(mns, opt, log10λ)
lsc = TransD_GP.testupdate(optlog10λ, log10λ)
ftest = TransD_GP.testupdate(opt, log10λ, mns)

f, ax = plt.subplots(2, 1, sharex=true)
ax[1].plot(xall', sqrt.(log10λ.fstar)')
ax[1].plot(xall', 10 .^lsc, "--k")
ax[1].plot(log10λ.xtrain[:,1:log10λ.n]', 10 .^(log10λ.ftrain[1:log10λ.n]), "+m")
ax[1].plot(log10λ.xtrain[:,log10λ.n]', 10^(log10λ.ftrain[log10λ.n]), ".r", markersize=10)
ax[1].set_title("λ")

ax[2].plot(xall', mns.fstar)
ax[2].plot(xall', ftest)
ax[2].plot(mns.xtrain[:,1:mns.n]', mns.ftrain[1:mns.n], "xr")
ax[2].plot(mns.xtrain[:,mns.n]', mns.ftrain[mns.n], ".m", markersize=10)
ax[2].set_title("function")
## update mns.fstar through death in mns
TransD_GP.death!(mns, opt)
lsc = TransD_GP.testupdate(optlog10λ, log10λ)
ftest = TransD_GP.testupdate(opt, log10λ, mns)

f, ax = plt.subplots(2, 1, sharex=true)
ax[1].plot(xall', sqrt.(log10λ.fstar)')
ax[1].plot(xall', 10 .^lsc, "--k")
ax[1].plot(log10λ.xtrain[:,1:log10λ.n]', 10 .^(log10λ.ftrain[1:log10λ.n]), "+m")
ax[1].plot(log10λ.xtrain[:,log10λ.n]', 10^(log10λ.ftrain[log10λ.n]), ".r", markersize=10)
ax[1].set_title("λ")

ax[2].plot(xall', mns.fstar)
ax[2].plot(xall', ftest)
ax[2].plot(mns.xtrain[:,1:mns.n]', mns.ftrain[1:mns.n], "xr")
ax[2].plot(mns.xtrain[:,mns.n]', mns.ftrain[mns.n], ".m", markersize=10)
ax[2].set_title("function")
## update λ and mns.fstar through death in λ
doall = false
lold = copy(log10λ.fstar)
TransD_GP.death!(log10λ, optlog10λ, mns, opt, doall=doall)
lsc = TransD_GP.testupdate(optlog10λ, log10λ)
ftest = TransD_GP.testupdate(opt, log10λ, mns)

f, ax = plt.subplots(2, 1, sharex=true)
ax[1].plot(xall', sqrt.(log10λ.fstar)')
ax[1].plot(xall', sqrt.(lold)')
ax[1].plot(xall', 10 .^lsc, "--k")
ax[1].plot(log10λ.xtrain[:,1:log10λ.n]', 10 .^(log10λ.ftrain[1:log10λ.n]), "+m")
ax[1].plot(log10λ.xtrain[:,log10λ.n]', 10^(log10λ.ftrain[log10λ.n]), ".r", markersize=10)
ax[1].set_title("λ")

ax[2].plot(xall', mns.fstar)
ax[2].plot(xall', ftest)
ax[2].plot(mns.xtrain[:,1:mns.n]', mns.ftrain[1:mns.n], "xr")
ax[2].set_title("function")

ax2 = ax[1].twinx()
ax2.plot(xall', log10.(abs.(lold' -10 .^2lsc)))
ax[1].grid()
ax[2].grid()
## update λ and mns.fstar through property change in λ
doall = false
lold = copy(log10λ.fstar)
TransD_GP.property_change!(log10λ, optlog10λ, mns, opt, doall=doall)
lsc = TransD_GP.testupdate(optlog10λ, log10λ)
ftest = TransD_GP.testupdate(opt, log10λ, mns)

f, ax = plt.subplots(2, 1, sharex=true)
ax[1].plot(xall', sqrt.(log10λ.fstar)')
ax[1].plot(xall', sqrt.(lold)')
ax[1].plot(xall', 10 .^lsc, "--k")
ax[1].plot(log10λ.xtrain[:,1:log10λ.n]', 10 .^(log10λ.ftrain[1:log10λ.n]), "+m")
ax[1].plot(log10λ.xtrain[:,log10λ.n]', 10^(log10λ.ftrain[log10λ.n]), ".r", markersize=10)
ax[1].set_title("λ")

ax[2].plot(xall', mns.fstar)
ax[2].plot(xall', ftest)
ax[2].plot(mns.xtrain[:,1:mns.n]', mns.ftrain[1:mns.n], "xr")
ax[2].set_title("function")

ax2 = ax[1].twinx()
ax2.plot(xall', log10.(abs.(lold' -10 .^2lsc)))
ax[1].grid()
ax[2].grid()
## update λ and mns.fstar through position change in λ
doall = false
lold = copy(log10λ.fstar)
TransD_GP.position_change!(log10λ, optlog10λ, mns, opt, doall=doall)
lsc = TransD_GP.testupdate(optlog10λ, log10λ)
ftest = TransD_GP.testupdate(opt, log10λ, mns)

f, ax = plt.subplots(2, 1, sharex=true)
ax[1].plot(xall', sqrt.(log10λ.fstar)')
ax[1].plot(xall', sqrt.(lold)')
ax[1].plot(xall', 10 .^lsc, "--k")
ax[1].plot(log10λ.xtrain[:,1:log10λ.n]', 10 .^(log10λ.ftrain[1:log10λ.n]), "+m")
ax[1].plot(log10λ.xtrain[:,log10λ.n]', 10^(log10λ.ftrain[log10λ.n]), ".r", markersize=10)
ax[1].set_title("λ")

ax[2].plot(xall', mns.fstar)
ax[2].plot(xall', ftest)
ax[2].plot(mns.xtrain[:,1:mns.n]', mns.ftrain[1:mns.n], "xr")
ax[2].set_title("function")

ax2 = ax[1].twinx()
ax2.plot(xall', log10.(abs.(lold' -10 .^2lsc)))
ax[1].grid()
ax[2].grid()
