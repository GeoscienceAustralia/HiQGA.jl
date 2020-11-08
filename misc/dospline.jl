using PyCall, PyPlot, Random
Random.seed!(5)
sp = pyimport("scipy.interpolate")
x = LinRange(0,1,50)
y = (0.9 .+ 0.1rand(length(x))).*sin.(2*pi*(x.-0.5))
t = collect(x[2:2:end-1])
s1 = sp.LSQUnivariateSpline(x, y, t)
x2 = LinRange(0, 1, 201) # new x-grid
y2 = s1(x2) # evaluate spline on that new grid
figure()
plot(x,y,label="original")
plot(x2,y2,label="interp", color="k")
knots = s1.get_knots()
c = s1.get_coeffs()
newknots(knots, k) = vcat(fill(knots[1],k),knots,fill(knots[end],k))
forscipyknots = newknots(knots, 3)
s2 = sp.BSpline(forscipyknots, c, 3)
y3 = s2(x2)
plot(x2,y3,"--r", label="reconstructed \nfrom knots and coeff")
legend()
##
knots = [0,.4,.4,.4,.4,.7,1]
c = [2,-5,5,2,-3,-1,2,1.5]
forscipyknots = newknots(knots, 3)
s2 = sp.BSpline(forscipyknots, c, 3)
figure()
plot(x2, s2(x2))
##
using Images
img = convert(Array{Float64, 2}, Gray.(load("func.png")))
