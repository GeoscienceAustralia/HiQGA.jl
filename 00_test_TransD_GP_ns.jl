Sys.iswindows() && (ENV["MPLBACKEND"]="qt4agg")
using PyPlot, Test, Random, Revise, Statistics, LinearAlgebra
any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
import GP, TransD_GP
## make options for the multichannel lengthscale GP
nminlog10λ, nmaxlog10λ = 2, 100
pnorm = 2.
Klog10λ = GP.Mat32()
λx,λy = 0.6,0.6
x = 0:(0.01λx):λx
y = 0:(0.01λy):2λy
λlog10λ = [0.1maximum(y), 0.1maximum(x)]
demean = true
sdev_poslog10λ = [0.01maximum(y), 0.01maximum(x)]
log10bounds = [-2 -0.69; -2 -0.69]
δlog10λ = 0.1
sdev_proplog10λ = [0.1, 0.1]
xall = zeros(2,length(x)*length(y))
for i in 1:size(xall,2)
    yid, xid = Tuple(CartesianIndices((length(y),length(x)))[i])
    xall[:,i] = [x[xid]; y[yid]]
end
xbounds = zeros(Float64,size(xall, 1), 2)
for dim in 1:size(xall, 1)
    xbounds[dim,:] = [minimum(xall[dim,:]), maximum(xall[dim,:])]
end
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
                        K = Klog10λ
                        )
@time  log10λ = TransD_GP.init(optlog10λ)
## make options for the nonstationary GP
nmin, nmax = 2, 400
fbounds = [-2. 2]
δ = 0.1
sdev_prop = [0.1]
sdev_pos = [0.05, 0.05]
K = GP.Mat32()
## Initialize model for the nonstationary GP
opt = TransD_GP.Options(nmin = nmin,
                        nmax = nmax,
                        xbounds = optlog10λ.xbounds,
                        fbounds = fbounds,
                        xall = optlog10λ.xall,
                        λ = λlog10λ,
                        δ = δ,
                        demean = demean,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        pnorm = pnorm,
                        quasimultid = false,
                        K = K
                        )
##
λ² = log10λ.fstar
@time m = TransD_GP.init(opt, log10λ)
idxs = TransD_GP.gettrainidx(opt.kdtree, m.xtrain, m.n)
ftest, = GP.GPfit(K, m.ftrain[:,1:m.n], m.xtrain[:,1:m.n],
    opt.xall, λ², λ²[:,idxs], δ, p=2, demean=demean, nogetvars=true)
@test norm(mean(ftest - m.fstar)) < 1e-12
