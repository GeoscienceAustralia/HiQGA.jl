using Documenter, Literate, Pkg
import HiQGA # 

# example_src_dir = joinpath(@__DIR__, "..", "examples", "1D", "stationary")
# example_src = joinpath(example_src_dir, "01_make_model.jl")
# example_file = joinpath(example_src_dir, "func2.txt")

# preprocess_1d(str) = replace(str, "func2.txt" => example_file)

outdir = joinpath(@__DIR__, "src", "generated")
# Literate.markdown(example_src, outdir; execute=true, documenter=true, preprocess=preprocess_1d)

makedocs(sitename="HiQGA - High Quality Geophysical Analysis",
         pages = [
                  "index.md",
#                 "Examples" => ["1D nonlinear regression" => ["generated/01_make_model.md", ],
#                               "SkyTEM" => []]
                 ]
        )
