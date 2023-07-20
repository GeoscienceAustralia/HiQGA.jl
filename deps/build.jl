ENV["PYTHON"]=""
using Pkg
Pkg.add("PyPlot")
Pkg.build("PyPlot")
error("Ooops")