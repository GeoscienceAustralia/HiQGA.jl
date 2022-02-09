using PyPlot
Tmax = 5
npoints = 10
g(n) = range(1/Tmax^n,1., length=npoints)
h(n) = range(0, n*log10(Tmax), length=npoints)
j(n) = range(0, 1-exp(-n*(Tmax-1)), length=npoints)

##
figure()
for n in range(0, 5, length=11)
    #plot(g(n).^(-1/n), g(n), label="T^$(-n)", ".-")
    plot(1 .- 1/n*log.(1 .- j(n)), j(n), "x-", label="1-exp(-$n(T -1))")
end
n = 1
plot(10 .^(1/n*h(n)), h(n), label="log T", "v-k")
xlabel("T")
ylabel("equally spaced f(T)")
legend(bbox_to_anchor=(1.04,1), borderaxespad=0)  
plt.tight_layout()      