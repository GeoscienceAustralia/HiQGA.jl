using transD_GP.GP, Random, PyPlot, Statistics, LinearAlgebra, PositiveFactorizations
##
Random.seed!(11)
f = 1.0/0.5
t = 0:(10*f)^-1:2/f
ntrain = 1
ntries = 20
λ = [0.1, 0.1]
δ = 0.01
δtry = 5δ
p = 2
demean = true
findmin = true
switchflag = findmin ? -1. : 1.
knownvalue = NaN
acqfun = GP.EI()

y = [100000*switchflag*(2x1+ sin.(2*pi*f*x1) + sinc.(2*f*(x1.-middle(x1))))*
                       (2x2 + sin.(2*pi*f*x2) + sinc.(2*f*(x2.-middle(x2)))) for x1 in t, x2 in t]
ynoisy = y + δ*randn(size(y))

X1, X2 = [x1 for x1 in t, x2 in t], [x2 for x1 in t, x2 in t]
T = [X1[:]';X2[:]']
tperm = randperm(length(y))
trainidx = tperm[1:ntrain]
ytrain = ynoisy[trainidx]
ttrain = T[:,trainidx]

vmin, vmax = extrema(ynoisy)
ytest, σ2, σ_prior  = GP.GPfit(GP.SqEuclidean(), ytrain', ttrain, T, λ.^2, δtry, p=p, demean=demean)
AF = GP.getAF(acqfun, ytrain, vec(ytest), diag(σ2), findmin=findmin, knownvalue=knownvalue)
f, ax = plt.subplots(3, ntries+1, figsize=(12,3), sharex=true, sharey=true)
ax[1].imshow(ynoisy, extent=[t[1], t[end], t[end], t[1]])
ax[1].plot(ttrain[2,:][:],ttrain[1,:][:], ".m",markersize=5)
ax[2].imshow(reshape(ytest, size(y)), extent=[t[1], t[end], t[end], t[1]], vmin=vmin, vmax=vmax)
ax[3].imshow(reshape(AF, size(y)), extent=[t[1], t[end], t[end], t[1]])
ax[1].set_title("0")
##
for i = 1:ntries
    global AF
    nextpos = argmax(AF)
    push!(ytrain, ynoisy[nextpos])
    ttrain = hcat(ttrain, T[:,nextpos])
    ytest, σ2, σ_prior  = GP.GPfit(GP.SqEuclidean(), ytrain', ttrain, T, λ.^2, δtry, p=p, demean=demean)
    AF = GP.getAF(acqfun, ytrain, vec(ytest), diag(σ2), findmin=findmin, knownvalue=knownvalue)
    ax[3i+1].imshow(y, extent=[t[1], t[end], t[end], t[1]])
    ax[3i+1].plot(ttrain[2,:][:],ttrain[1,:][:], ".m",markersize=5)
    ax[3i+2].imshow(reshape(ytest, size(y)), extent=[t[1], t[end], t[end], t[1]], vmin=vmin, vmax=vmax)
    ax[3i+3].imshow(reshape(AF, size(y)), extent=[t[1], t[end], t[end], t[1]])
    ax[3i+1].set_title("$i")
end
plt.tight_layout()