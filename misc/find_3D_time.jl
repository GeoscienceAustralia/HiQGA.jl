Sys.iswindows() && (ENV["MPLBACKEND"]="qt4agg")
using PyPlot, Test, Random, Revise, Statistics, LinearAlgebra
srcdir = dirname(pwd())*"/src"
any(srcdir .== LOAD_PATH) || push!(LOAD_PATH, srcdir)
import GP, TransD_GP
## make options for the multichannel lengthscale GP
nminlog10λ, nmaxlog10λ = 2, 200
pnorm = 2.
Klog10λ = GP.Mat32()
λx,λy,λz = 1, 1, 1
x = 0:(0.013λx):λx
y = 0:(0.013λy):λy
z = 0:(0.013λy):λy
λlog10λ = [0.05maximum(x), 0.05maximum(y), 0.05maximum(z)]
demean = false
sdev_poslog10λ = [0.01maximum(y), 0.01maximum(x), 0.01maximum(z)]
log10bounds = [-1 -0.69; -1 -0.69; -1 0.69]
δlog10λ = 0.1
sdev_proplog10λ = [0.1, 0.1, 0.1]
xall = zeros(3,length(x)*length(y)*length(z))
for i in 1:size(xall,2)
    xid, yid, zid = Tuple(CartesianIndices((length(x),length(y), length(z)))[i])
    xall[:,i] = [x[xid]; y[yid]; z[zid]]
end
xbounds = zeros(Float64,size(xall, 1), 2)
for dim in 1:size(xall, 1)
    xbounds[dim,:] = [minimum(xall[dim,:]), maximum(xall[dim,:])]
end
## Initialize a model using these options
Random.seed!(12)
optlog10λ = TransD_GP.OptionsStat(nmin = nminlog10λ,
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
                        timesλ = 3.6
                        )
@time  log10λ = TransD_GP.init(optlog10λ)
## make options for the nonstationary GP
nmin, nmax = 2, 200
fbounds = [-2. 2]
δ = 0.1
sdev_prop = [0.1]
sdev_pos = [0.05, 0.05, 0.05]
K = GP.Mat32()
demean_ns = true
## Initialize model for the nonstationary GP
Random.seed!(13)
opt = TransD_GP.OptionsNonstat(optlog10λ,
                        nmin = nmin,
                        nmax = nmax,
                        fbounds = fbounds,
                        δ = δ,
                        demean = demean_ns,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        pnorm = pnorm,
                        K = K
                        )
@time m = TransD_GP.init(opt, log10λ)
## timing for s birth in s model for ns mode
log10λ = TransD_GP.init(optlog10λ)
m = TransD_GP.init(opt, log10λ)
for i = 1:98
    TransD_GP.birth!(log10λ, optlog10λ, m, opt)
end
@time for i = 98
    TransD_GP.birth!(m, opt, log10λ)
end
NTIMES = 1
ntimes = 150
T = zeros(NTIMES)
for I = 1:NTIMES
    T[I] = time()
    for i = 1:ntimes
        TransD_GP.death!(log10λ, optlog10λ, m, opt)
        TransD_GP.birth!(log10λ, optlog10λ, m, opt)
    end
    T[I] = (time() - T[I])/2ntimes
end
@info "time for $ntimes birth/death is $(mean(T)) +- $(std(T)/sqrt(NTIMES))"
