# # Example - Bayesian nonlinear regression
# ## Setup
## let's try gradient descent, all model values are in log10 conductivity
ρstart, ρ0 = map(y->zeros(length(x)), 1:2)
ρstart .= 0.
ρ0 .= 0
regtype = :R1
lo = -5
hi = 5
## one step
A = line.J
d = line.d[line.select]
δ² = 1e-6
R = transD_GP.makereg(:R1, line)
m_mle = (A'A + δ²*R'R)\(A'd)

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
lalpha = 0.1
idxlast=length(m)
for (i, mi) in enumerate(m)
    if isempty(λ²[i]) 
        global idxlast = i-1
        break
    end    
    f = figure(figsize=(7,6))
    s1 = subplot(131)
    for ii in 1:length(χ²[i])
        s1.plot(mi[ii], x) 
    end
    s1.plot(mi[idx[i]], x, color="r", linewidth=2)
    s1.plot(y, x, color="k", linewidth=2, linestyle="--")
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
for (i, mi) in enumerate(m)
    plot(mi[idx[i]], x, alpha=alpha)
    i == idxlast && break
end
plot(m[idxlast][idx[idxlast]], x, color="r", linewidth=2)
plot(y, x, color="k", linewidth=2, linestyle="--")
gca().invert_yaxis()
gca().invert_xaxis()
ylabel("depth m")
xlabel("Log₁₀ ρ")
plt.tight_layout()
ax = gca()
## Compare with posterior model covariance
using LinearAlgebra
F = line
R = transD_GP.makereg(regtype, F)
JtW, Wr = F.J'*F.W, F.W*F.res
H = JtW*(JtW)' + λ²[idxlast][idx[idxlast]][1]*R'R
Cpost = inv(Hermitian(Matrix(H)))
figure()
zplot = [z; z[end] + 5]
pcolormesh(zplot, zplot, Cpost)
xlabel("Depth m")
ylabel("Depth m")
colorbar()
sd = sqrt.(diag(Cpost))
ax.fill_betweenx(zplot[1:end-1]+diff(zplot)/2, m[idxlast][idx[idxlast]] -sd, m[idxlast][idx[idxlast]] +sd, alpha=0.2)
ax.set_xlim(lo, hi)