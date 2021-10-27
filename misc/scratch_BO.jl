using transD_GP.GP, Random, PyPlot, Statistics, LinearAlgebra, PositiveFactorizations
##
Random.seed!(12)
f = 1.0/0.5
t = 0:(100*f)^-1:2/f
ntrain = 1
ntries = 4
λ = [0.1]
δ = 0.01
δtry = 5δ
p = 2
demean = true
findmin = true
switchflag = findmin ? -1. : 1.
knownvalue = NaN
acqfun = GP.EI()

y = 100000*switchflag*(2t+ sin.(2*pi*f*t) + sinc.(2*f*(t.-middle(t))))
ynoisy = y + δ*randn(length(t))

tperm = randperm(length(t))
trainidx = tperm[1:ntrain]
ytrain = ynoisy[trainidx]
ttrain = t[trainidx]

ytest, σ2, σ_prior  = GP.GPfit(GP.SqEuclidean(), ytrain', ttrain', t', λ.^2, δtry, p=p, demean=demean)
AF = GP.getAF(acqfun, ytrain, vec(ytest), diag(σ2), findmin=findmin, knownvalue=knownvalue)
f, ax = plt.subplots(2, ntries+1, figsize=(12,5), sharex="col", sharey="row")
ax[1].plot(t,y, "-k")
ax[1].plot(t,ytest)
ax[1].plot(ttrain,ytrain, ".m",markersize=20)
ax[1].plot(t,ytest+2sqrt.(diag(σ2)),"--r",alpha=0.2)
ax[1].plot(t,ytest-2sqrt.(diag(σ2)),"--b",alpha=0.2)
ax[1].grid()
ax[2].plot(t, AF)
ax[2].grid()
ax[1].set_title("step 0")
##
for i = 1:ntries
    global AF
    nextpos = argmax(AF)
    push!(ytrain, ynoisy[nextpos])
    push!(ttrain, t[nextpos])
    ytest, σ2, σ_prior  = GP.GPfit(GP.SqEuclidean(), ytrain', ttrain', t', λ.^2, δtry, p=p, demean=demean)
    AF = GP.getAF(acqfun, ytrain, vec(ytest), diag(σ2), findmin=findmin, knownvalue=knownvalue)
    ax[2i+1].plot(t,y, "-k")
    ax[2i+1].plot(ttrain,ytrain, ".m",markersize=20)
    ax[2i+1].plot(t,ytest)
    ax[2i+1].plot(t,ytest+2sqrt.(diag(σ2)),"--r",alpha=0.2)
    ax[2i+1].plot(t,ytest-2sqrt.(diag(σ2)),"--b",alpha=0.2)
    ax[2i+1].grid()
    ax[2(i+1)].plot(t, AF)
    ax[2(i+1)].grid()
    ax[2i+1].set_title("step $i")
end
# savefig("fig.png",dpi=300)

##
figure()
s1=subplot(211, aspect="auto")
imshow(σ_prior, extent=[t[1],t[end],t[end],t[1]])
s2=subplot(212, sharex=s1, sharey=s1)
imshow(σ2, extent=[t[1],t[end],t[end],t[1]])

figure()
plot(t,y, "-k",label="true")
plot(ttrain,ytrain, "+m", label="train",markersize=10)
nsamp = length(t)
nreal=50
realizations = randn(nreal,nsamp)*cholesky(Positive, var(ytest)*σ2 + 1e-12*Matrix(I,nsamp,nsamp)).U
realizations = ytest .+ realizations'
plot(t,realizations,alpha=0.2)
