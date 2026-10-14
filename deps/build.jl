import Pkg, Conda
@info "building!"

# Using conda-forge; the Anaconda `defaults` channel has licensing/solver problems.
Conda.add_channel("conda-forge")

# PyCall is not GC-safe on Python 3.12+ (segfaults in Py_Dealloc). Pin to 3.11.
Conda.add("python=3.11")

Conda.add("matplotlib")

# "" tells PyCall to use Conda.jl's root environment and records PyCall.conda = true,
# which lets PyPlot auto-install matplotlib if it's ever missing.
ENV["PYTHON"] = ""
Pkg.build("PyCall")

@info "built PyCall!"