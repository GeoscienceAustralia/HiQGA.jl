import Pkg
Pkg.add("Conda")
using Conda
Conda.add("matplotlib")
ENV["PYTHON"] = joinpath(Conda.ROOTENV, "bin", "python")
Pkg.add("PyCall"); Pkg.build("PyCall")
Pkg.add("PyPlot")';