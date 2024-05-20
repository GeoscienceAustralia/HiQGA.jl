import Pkg, Conda
@info "building!"
Conda.pip_interop(true)
Conda.pip("install", "matplotlib==3.7.3")
Conda.add("matplotlib")
ENV["PYTHON"] = joinpath(Conda.ROOTENV, "bin", "python")
Pkg.build("PyCall")
@info "built PyCall!"
