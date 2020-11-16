# transD_GP.jl

![CI status](https://github.com/a2ray/transD_GP/workflows/CI/badge.svg)

This package implements the trans-dimensional Gaussian process algorithm as published in [*Bayesian geophysical inversion with trans-dimensional Gaussian process machine learning*, A. Ray and D. Myer, Geophysical Journal International **217(3)**, 2019](https://doi.org/10.1093/gji/ggz111).

## Installation
To install, use Julia's `Pkg` REPL:
```
pkg> add https://github.com/a2ray/transD_GP.git
```

## Usage
Examples of how to use the package can be found in the `examples` directory. The Markov Chain Monte Carlo sampler is configured to support parallel tempering on multiple CPUs - some of the examples accomplish this with Julia's built-in multiprocessing, and others use MPI in order to support inversions on HPC clusters that don't work with Julia's default SSH-based multiprocessing. The MPI examples require [MPI.jl](https://github.com/JuliaParallel/MPI.jl) and [MPIClusterManagers.jl](https://github.com/JuliaParallel/MPIClusterManagers.jl/), which are not installed as dependencies for this package, so you will need to ensure they are installed and configured correctly to run these examples. Please note that MPIClusterManagers.jl has issues with Julia <1.4.2, so please ensure you are using an up-to-date Julia version.