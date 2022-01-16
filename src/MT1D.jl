module MT1D
using PyPlot
const μ = 4*pi*1e-7

function Z_f(freqs, ρ, h)
    n = length(ρ)
    @assert length(h) == n - 1
    [Z_f(f, ρ, h, n) for f in freqs]   
end    

function Z_f(f, ρ, h, n)
    # from Kaufman & Keller 1981 
    ω = 2pi*f
    k = sqrt(im*ω*μ/ρ[end])
    Z = ω*μ/k
    for i = n-1:-1:1
        k = sqrt(im*ω*μ/ρ[i])
        ωμ_over_k = ω*μ/k
        Z = ωμ_over_k*coth(-im*k*h[i] + acoth(Z/ωμ_over_k))
    end
    Z    
end

ρapp(freqs, Z) = abs2.(Z)./(2pi*freqs*μ)

phase(Z) = rad2deg.(angle.(Z))

function twolayer_ex(h1, ρ1, ρ2::AbstractArray; ntimesperdecade=10)
    X = 10 .^(-1:1/ntimesperdecade:5) # λ₁/h₁
    T = ((h1*X/(2pi*503)).^2)/ρ1
    Z = [twolayer_ex(h1, ρ1, rho2, T) for rho2 in ρ2]
    fig, ax = plt.subplots(2,1, sharex=true, figsize=(6,10))
    for ZZ in Z
        ax[1].loglog(X, ρapp(1 ./T, ZZ)/ρ1)
        ax[2].semilogx(X, phase(ZZ))
    end
    ax[1].set_ylim(1e-3, 1e4)
    ax[2].set_ylim(-95, 10)
    ax[1].grid()
    ax[2].grid()
    ax[2].set_xlabel(L"\lambda_1/h_1")
    ax[1].set_ylabel(L"\rho_a/\rho_1")
    ax[2].set_ylabel("Phase "*L"^\circ")
    fig.suptitle("Kaufman & Keller 1981 Fig 3.1")
    fig.tight_layout()  
    nothing
end

function twolayer_ex(h1, ρ1, ρ2::S, T) where S<:Number
    ρ = [ρ1, ρ2]
    Z_f(1 ./T, ρ, [h1])
end    

function plotcurve(T, Z; showfreq=false, gridalpha=0.5)
    fig, ax = plt.subplots(1, 2, sharex=true)
    plotcurve(T, Z, fig, showfreq=showfreq, gridalpha=gridalpha)
    fig
end    

function plotcurve(T, Z, fig; showfreq=false, iaxis=1, gridalpha=0.5) 
    f = 1 ./T
    ρₐ = ρapp(f, Z)
    ϕ  = phase(Z)
    if showfreq 
        abcissa = f
        xlabel = "Frequency Hz"
    else
        abcissa = T
        xlabel = "Time s"
    end    
    ax = fig.axes
    ax[iaxis].loglog(abcissa, ρₐ)
    ax[iaxis].set_xlabel(xlabel)
    ax[iaxis].set_ylabel(L"\rho_{app}"*" ohm-m")
    ax[iaxis].grid(b=true, which="both", alpha=gridalpha)
    ax[iaxis+1].semilogx(abcissa, ϕ)
    ax[iaxis+1].set_xlabel(xlabel)
    ax[iaxis+1].set_ylabel("Phase "*L"^\circ")
    ax[iaxis+1].grid(b=true, which="both", alpha=gridalpha)
    fig.tight_layout()
end    

function plotmodelcurve(T, ρ, z; showfreq=false, figsize=(10,4), gridalpha=0.5)
    fig = figure(figsize=(figsize))
    s1 = subplot(131)
    s2 = subplot(132)
    s3 = subplot(133, sharex=s2)
    plotmodelcurve(T, ρ, z, fig, showfreq=showfreq, gridalpha=gridalpha)
    fig
end

function plotmodelcurve(T, ρ, z, fig; showfreq=false, gridalpha=0.5)
    f = 1 ./T
    h = diff(z)
    Z = Z_f(f, ρ, h)
    ax = fig.axes
    zlast = diff(z)[end] +z[end]
    ax[1].step([ρ; ρ[end]], [z; zlast])
    ax[1].set_xlabel("ρ ohm-m")
    ax[1].set_ylabel("Depth m")
    y1, y2 = ax[1].get_ylim()
    y1 < y2 && ax[1].invert_yaxis()
    ax[1].set_xscale("log")
    ax[1].grid(b=true, which="both", alpha=gridalpha)
    plotcurve(T, Z, fig, showfreq=showfreq, iaxis=2, gridalpha=gridalpha) 
end

end