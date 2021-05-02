using PyPlot, Random, Revise, Statistics, LinearAlgebra,
      Distributed, DelimitedFiles, transD_GP
## 1D functions
Random.seed!(200)
A = readdlm("hp.csv")[1:1:end,:]
idx = sortperm(A[:,1])
t = A[idx,1]
ly = log10.(A[idx,2])
function bindata(x, y; ninbinmin=10, binwidth=10)
    Y, σy, X = zeros(0), zeros(0), zeros(0)
    @assert length(x) == length(y)
    n = length(x)
    pointcount = 1
    while pointcount<n
         sum_x  = 0.
         sum_y  = 0.
         sum_y² = 0.
         ninbin = 0
         δx     = 0.
         agestart = x[pointcount]
         while true
             sum_y  += y[pointcount]
             sum_y² += y[pointcount]^2
             sum_x  += x[pointcount]
             δx = x[pointcount] - agestart
             ninbin += 1
             pointcount += 1
             ((ninbin >= ninbinmin) & (δx >= binwidth)) && break
             pointcount == n && break
         end
         @info ninbin
         push!(Y,  sum_y/ninbin)
         # push!(??y, sqrt((sum_y�?/ninbin - (sum_y/ninbin)^2)/ninbin))
         push!(σy, sqrt((sum_y²/ninbin - (sum_y/ninbin)^2)))
         if ninbin == 1 && pointcount == n
             σy[end] = σy[end-1]
         end
         # push!(X, sum_x/ninbin)
         push!(X, agestart+δx/2)
    end
    X, Y, σy
end
##
X, Y, σy = bindata(t, ly, ninbinmin=10, binwidth=10)
figure()
errorbar(X, Y, yerr = σy, marker=".", elinewidth=1, capsize=3)
xlabel("time in My")
ylabel("f(t)")
plt.tight_layout()
line = transD_GP.Line(Y;useML=false, σ=σy)
