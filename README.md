# transD_GP.jl

![CI status](https://github.com/a2ray/transD_GP/workflows/CI/badge.svg)

This package implements both the nested (2-layer) and vanilla trans-dimensional Gaussian process algorithm as published in 
- [*Bayesian inversion using nested trans-dimensional Gaussian processes*, A. Ray, Geophysical Journal International, **226(1)**, 2021](https://doi.org/10.1093/gji/ggab114).
- [*Bayesian geophysical inversion with trans-dimensional Gaussian process machine learning*, A. Ray and D. Myer, Geophysical Journal International **217(3)**, 2019](https://doi.org/10.1093/gji/ggz111).

Readily usable geophysical forward operators are to do with CSEM and AEM physics (references underneath), **for which the time domain AEM codes are fairly production-ready**.

## Installation
To install, in a perfect world we'd use Julia's `Pkg` REPL by hitting `]` to enter `pkg>` mode. Then enter the following, at the `pkg>` prompt:
```
pkg> add https://github.com/a2ray/transD_GP.git
```
***But we can't do this anymore with GitHub requiring the use of ssh-keys as password authentication has been/soon will be outdated.*** Instead we now follow the instructions here to [create your own ssh-keys](https://docs.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent) and then [adding your ssh-keys to Github](https://docs.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account). You don't need to do this if you're already using ssh-keys with GitHub, and can simply follow on from here. You should create the keys on the computer you'd like to use `transD_GP` from, most likely your laptop or on your friendly local cluster. Do the following at your shell prompt:
```
cd ~/.ssh
eval $(ssh-agent)
ssh-add id_ed25519
ssh-add -l
julia
```
at the REPL you go to `Pkg` mode by hitting `]` and then at the `pkg>` prompt:
```
pkg> add git@github.com:a2ray/transD_GP.git
```
If you get an error message saying the public key wasn't found, well, you'll have to add the public key `id_ed25519.pub` residing on GitHub into the abovementioned `.ssh` directory as well. Don't ever share your private key!!! All this pain will disappear and the perfect world instructions will apply if we make the repository public. We're almost, but not quite ready to do that yet.

## Usage
Examples of how to use the package can be found in the `examples` directory. The Markov Chain Monte Carlo sampler is configured to support parallel tempering on multiple CPUs - some of the examples accomplish this with Julia's built-in multiprocessing, and others use MPI in order to support inversions on HPC clusters that don't work with Julia's default SSH-based multiprocessing. The MPI examples require [MPI.jl](https://github.com/JuliaParallel/MPI.jl) and [MPIClusterManagers.jl](https://github.com/JuliaParallel/MPIClusterManagers.jl/), which are not installed as dependencies for this package, so you will need to ensure they are installed and configured correctly to run these examples. Please note that MPIClusterManagers.jl has issues with Julia <1.4.2, so please ensure you are using an up-to-date Julia version. 

Some example scripts have as a dependency [Revise.jl](https://github.com/timholy/Revise.jl) as we're still actively [developing this package](https://julialang.github.io/Pkg.jl/v1.5/getting-started/), so you may need to install [Revise](https://github.com/timholy/Revise.jl) if not already installed. All Julia users should be developing with [Revise](https://github.com/timholy/Revise.jl)! After installation, to run the examples, simply clone the package separately (or download as a ZIP), navigate to the `examples` folder and run the scripts in their numerical order.

## Updating
You'll have to run `ssh-agent` as shown above for the installation, but instead of adding the package again in `pkg>`, you would simply do at the `pkg>` prompt
```
pkg> update transD_GP
```

## Developing transD_GP or modifying it for your own special forward physics
```
pkg>dev git@github.com:a2ray/transD_GP.git
```
All the `ssh-keys` jazz above still holds. It will download to your `JULIA_PKG_DEVDIR`. [Here's a gist](https://gist.github.com/a2ray/8c2c55c25fee6647501b403886bbe64d) on adding your own module if you want to modify the source code. Alternatively, if you only want to use the sampling methods in `transD_GP` without contributing to the source (boo! j/k) [here's another gist](https://gist.github.com/a2ray/92a8c14483c21dda6ddf56685b95fbb8) which is more appropriate.

### References for AEM and CSEM physics 

- [Blatter, D., Key, K., Ray, A., Foley, N., Tulaczyk, S., & Auken, E. (2018). Trans-dimensional Bayesian inversion of airborne transient EM data from Taylor Glacier, Antarctica. Geophysical Journal International, 214(3)](https://doi.org/10.1093/gji/ggy255)

- [Ray, A., & Key, K. (2012). Bayesian inversion of marine CSEM data with a trans-dimensional self parametrizing algorithm. Geophysical Journal International, 191(3), 1135-1151.](https://doi.org/10.1111/j.1365-246X.2012.05677.x)
