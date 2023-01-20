# # Example - Bayesian nonlinear regression
# ## Setup
using PyPlot, Random, LinearAlgebra, SparseArrays
using HiQGA.transD_GP
## true model
zstart = 0
extendfrac, dz = 1.03, 1.5
zall, znall, zboundaries, = transD_GP.setupz(zstart, extendfrac, dz=dz, n=65)
z, ρ, nfixed = transD_GP.makezρ(zboundaries; zfixed=[], ρfixed=[])
# fill in detail in ohm-m
ρ[(z.>=zstart) .& (z.<50)] .= 20.
ρ[(z.>=50) .&(z.<80)] .= 1
ρ[(z.>=80) .&(z.<100)] .= 20
ρ[(z.>=100) .&(z.<200)] .= 50
ρ[(z.>=200) .&(z.<250)] .= 80
ρ[(z.>=250)]            .= 150
# add jitter to model in log10 domain
Random.seed!(11)
ρ = 10 .^(0.1*randn(length(ρ)) + log10.(ρ))
## randomly sample this and add Gaussian noise of sdev sigma (in log10)
σ, fracdrop = 0.2, 0.4 
Random.seed!(9)
ρnoisy = 10 .^(σ*randn(length(ρ)) + log10.(ρ))
ndrop = round(Int, fracdrop*length(ρ))
linidx = randperm(length(ρ))[1:ndrop]
ρnoisy[linidx] .= NaN 
# make the operator of type Line
line = transD_GP.Line(log10.(ρnoisy), σ=σ, calcjacobian=true)
## Plot
f = figure(figsize=(16,8))
# step(1:length(z), log10.(ρ))
errorbar(znall, log10.(ρnoisy), yerr=σ, ls="none", marker=".", capsize=8, markersize=32)
xlabel("x")
ylabel("f(x)")
grid(axis="x")
gca().invert_yaxis()
transD_GP.nicenup(gcf(), fsize=18)
# plot system Matrix
figure()
imshow(line.J, cmap="hot", extent=[0.5,65.5,39.5, 0.5])
transD_GP.nicenup(gcf(), fsize=18)
## do the inversion
## let's try gradient descent, all model values are in log10 conductivity
ρstart, ρ0 = map(x->zeros(length(zall)), 1:2)
ρstart .= 0.5
ρ0 .= 0
regtype = :R1
lo = -1
hi = 3
## one step
A = line.J
d = line.d[line.select]
δ² = 1e12
R = transD_GP.makereg(regtype, line)
m_mle = (A'A + δ²*R'R)\(A'd)
m_ridge = (A'A + δ²*I)\(A'd)
## do it
m, χ², λ², idx = transD_GP.gradientinv(ρstart, ρ0, line, nstepsmax=10, 
                            regularizeupdate=false, 
                            λ²min=-2, 
                            λ²max=7,
                            β² = 0.0,
                            ntries=15,
                            lo = lo,
                            hi = hi,
                            knownvalue=NaN, 
                            breakonknown = true,
                            firstvalue=:last, 
                            regtype = regtype);
## debug plots: all in each
alpha = 0.1
for (i, mi) in enumerate(m)
    f = figure(figsize=(7,6))
    s1 = subplot(131)
    for ii in 1:length(χ²[i])
        s1.step(mi[ii], z) 
    end
    s1.step(mi[idx[i]], z, color="r", linewidth=2)
    s1.step(log10.(ρ), z, color="k", linewidth=2, linestyle="--")
    s1.set_ylabel("depth m")
    s1.set_xlabel("Log₁₀ ρ")
    s1.invert_yaxis()
    s1.invert_xaxis()
    s2 = subplot(232)
    hps = hcat(λ²[i]...)'
    sortedλidx = sortperm(hps[:,1])
    sortedαidx = sortperm(hps[:,2])
    s2.loglog(hps[sortedλidx,1], χ²[i][sortedλidx], ".-", markersize=5)
    s2.plot(hps[idx[i],1], χ²[i][idx[i]], "o" )
    s2.set_ylabel("χ²")
    s3 = subplot(236)
    s3.loglog(χ²[i][sortedαidx], hps[sortedαidx,2], ".-", markersize=5)
    s3.plot(χ²[i][idx[i]], hps[idx[i],2], "o" )
    s3.set_xlabel("χ²")
    s4 = subplot(235, sharex=s2, sharey=s3)
    s4.scatter(hps[:,1], hps[:,2], c=log10.(χ²[i]), cmap="magma")
    s4.plot(hps[idx[i],1], hps[idx[i],2], "x", markersize=10)
    s4.set_xlabel("λ²"); s4.set_ylabel("step length α")
    plt.tight_layout()
end
## debug plots best in each
figure(figsize=(3,6))
alpha = 0.1
idxlast=length(m)
step(m[idxlast][idx[idxlast]], z, color="r", linewidth=2)
step(log10.(ρ), z, color="k", linewidth=2, linestyle="--")
gca().invert_yaxis()
gca().invert_xaxis()
ylabel("depth m")
xlabel("Log₁₀ ρ")
plt.tight_layout()
## compare MLE, truth and Occam
ocm = m[idxlast][idx[idxlast]]    
figure(figsize=(12,4))
step(1:length(z), m_mle, color="orange", 
    label=L"x_{smooth}^{MLE}"*" distance from truth = $(round((norm(m_mle-log10.(ρ))),sigdigits=2))")
step(1:length(z), ocm, color="r", linewidth=2, 
    label=L"x_{occam}"*" distance from truth = $(round((norm(ocm-log10.(ρ))), sigdigits=2))")
step(1:length(z), log10.(ρ), color="k", linewidth=2, label="truth")
errorbar(znall, log10.(ρnoisy), yerr=σ, ls="none", marker=".", capsize=8, markersize=16, alpha=0.25)
gca().invert_yaxis()
legend()
xlabel("x")
ylabel("f(x)")
transD_GP.nicenup(gcf(), fsize=16)
## try coordinate descent
using Wavelets, LinearMaps
wt = wavelet(WT.haar)
# wt = wavelet(WT.cdf97, WT.Lifting)
ocm = m[idxlast][idx[idxlast]][1:end-1] # must be a power of 2 to be useful ...
struct Fop
    wt
    n
end
struct Ftop
    wt
    n
end
Fop(;x=[1.], wt = wavelet(WT.haar)) = Fop(wt, length(x))
Ftop(foo::Fop) = Ftop(foo.wt, foo.n)
(foo::Fop)(x) = dwt(x, foo.wt)
(foo::Ftop)(x) = idwt(x, foo.wt)
F = Fop(;x=ocm, wt)
Ft = Ftop(F)
L = LinearMap(F, Ft, F.n)
A = transD_GP.LineRegression.getA(line.d)[:,1:end-1]*L'
# x = zeros(size(A,2))
x = L(ocm)
transD_GP.coordinatedesc(Matrix(A), x, d, 10 .^range(5, -5, 500), line.σ[line.select])
## fooling with trend filtering
A = transD_GP.LineRegression.getA(line.d)
C = transD_GP.makeinverseR1(size(A,2); η=1e0)
zz = zeros(size(A,2))
transD_GP.coordinatedesc(A*C, zz, d, 10 .^range(6, -6, 500), line.σ[line.select])
xx = C*zz
## plot
figure()
plot(z, log10.(ρ), label="true")
plot(z[1:end-1], L'(x), label="wt")
plot(z[1:end-1], ocm, label="occam")
plot(z, xx, label="TV")
legend()