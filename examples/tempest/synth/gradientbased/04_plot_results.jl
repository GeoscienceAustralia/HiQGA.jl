## one plot for each outer iteration
using Printf
alpha = 0.3
ndata = length(aem.res)
for (i, mi) in enumerate(m)
    f, s = plt.subplots(1, 2, gridspec_kw=Dict("width_ratios" => [1,1.5]),
        figsize=(8,5))
    for ii in 1:length(χ²[i])
        s[1].step(-mi[ii], aem.z[2:end], alpha=alpha) 
    end
    s[1].step(-mi[idx[i]], aem.z[2:end], color="r", linewidth=2)
    s[1].step(log10.(ρ[2:end]), z[2:end], color="k", linewidth=2, linestyle="--")
    s[1].set_ylabel("depth m")
    s[1].set_xlabel("Log₁₀ ρ")
    s[1].invert_yaxis()
    s[1].invert_xaxis()
    hps = hcat(λ²[i]...)'
    sortedλidx = sortperm(hps[:,1])
    sortedαidx = sortperm(hps[:,2])
    s[2].loglog(hps[sortedλidx,1], χ²[i][sortedλidx], ".-", markersize=5)
    s[2].plot(hps[idx[i],1], χ²[i][idx[i]], "o" )
    s[2].text(.95hps[idx[i],1], 1.05χ²[i][idx[i]], "ϕd = $(@sprintf("%.2f", χ²[i][idx[i]]/ndata)), χ²=$(@sprintf("%.2f", χ²[i][idx[i]]))" )
    s[2].plot(10. .^[λ²min, λ²max], [ndata, ndata], "--k")
    s[2].set_ylabel("χ²")
    s[2].set_xlabel("λ²")
    plt.suptitle("Iteration $i, α=$(hps[idx[i],2])")
    transD_GP.nicenup(gcf(), fsize=12)
end
## plots for best in each
_, s = plt.subplots(1, 2, gridspec_kw=Dict("width_ratios" => [1,1.5]),
        figsize=(8,5))
alpha = 0.3
χ²best = zeros(0)
idxlast = length(m)
ndata = length(aem.res)
for (i, mi) in enumerate(m)
    s[1].step(-mi[idx[i]], aem.z[2:end], alpha=alpha)
    push!(χ²best, minimum(χ²[i]))
    i == idxlast && break
end
s[1].step(-m[idxlast][idx[idxlast]], aem.z[2:end], color="r", linewidth=2)
s[1].step(log10.(ρ[2:end]), z[2:end], color="k", linewidth=2, linestyle="--")
s[1].invert_yaxis()
s[1].invert_xaxis()
s[1].set_ylabel("depth m")
s[1].set_xlabel("Log₁₀ ρ")
s[2].plot(1:idxlast, χ²best/ndata, ".-", ms=10)
s[2].plot([1, idxlast], [1, 1], "--k")
s[2].set_xlabel("Outer iteration #")
s[2].set_ylabel(L"ϕ_d"*" best")
gcf().tight_layout()