This folder contains various utilities and scripts as well as the package RDP.jl (Ramer-Douglas-Peucker) for GA's Earthsci and viewing portals, as well as reprojecting, generating shape files, making SEG-Y etc. Some of these depend on GDAL and other external libraries which would make HiQGA too big if registered all together. To install RDP and use these scripts, please first install HiQGA, find your HiQGA installation path, and develop RDP.jl like so:
```julia
using HiQGA, Pkg
pathdir = dirname(dirname(pathof(HiQGA))) # to get to root of HiQGA install directory
Pkg.develop(;path=joinpath(pathdir, "zz_portalcurtains") # which has its own src directory
```
Please note, RDP.jl is not a registered julia package, so anytime you update HiQGA, if any changes have been made to registered HiQGA and the included RDP src directory, changes will propagate here.