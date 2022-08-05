using StatsBase, Distributions
Random.seed!(13)
f, ax = plt.subplots(18, 3, gridspec_kw=Dict("width_ratios" => [1,5,1]),
        figsize=(15,12), sharey=true)
n = 15
s = zeros(size(ax, 1))
for i in 1:size(ax, 1)
    global s
    x = randn(n)
    h = fit(Histogram, x)
    height = diff(h.edges[1])
    ax[i,1].barh(h.edges[1][1:end-1], h.weights, align="edge", height=height)
    ax[i,2].sharey(ax[1])
    ax[i,2].plot(1:length(x),x,".-r", markersize=7)
    ax[i,2].plot(1:length(x),zeros(length(x)),"--k")
    ax[i,1].set_xticklabels([])
    ax[i,3].tick_params(labelbottom=false, labelleft=false, labelright=false, labeltop=false,
                    bottom=false, left=false, right=false, top=false)
    s[i] = round(x'x; sigdigits=3)
    ax[i,3].annotate("$(s[i])", xy=[0.5;0.5], 
            xycoords="axes fraction", ha="center", fontsize=12)
    i != size(ax, 1) && ax[i,2].set_xticklabels([])
    ax[i,2].set_xlim(1,length(x))
    ax[i,1].set_ylim(-3,3)
    if i==1
        ax[i,1].set_title("p(x)")
        ax[i,2].set_title("x = $n samples from N(0,1)")
        ax[i,3].set_title(L"\chi^2 = \sum_i x_i^2")
    end        
end
transD_GP.nicenup(f, fsize=12)
##
f = figure(figsize=(10,5))
h = fit(Histogram, s)
width = diff(h.edges[1])
w = (h.weights/sum(h.weights))./(width)
bar(h.edges[1][1:end-1], w, width, align="edge")
x = 1:n+40
p = pdf.(Chisq(n), x)
plot(x, p, "--k")
xlabel("χ²")
ylabel("p(χ²)")
axn = gca().twiny()
ax.get_shared_x_axes().join(ax, axn)
axn.plot(x/n, p, "--k")
axn.xaxis.label.set_color("red")
axn.tick_params(axis="x", colors="red")
axn.set_xlabel("χ²/$n")
transD_GP.nicenup(f, fsize=14)