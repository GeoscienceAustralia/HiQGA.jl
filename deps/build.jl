import Pkg, Conda
@info "building!"
Conda.add("matplotlib")
ENV["PYTHON"] = joinpath(Conda.ROOTENV, "bin", "python")
Pkg.build("PyCall")
@info "built PyCall!"
