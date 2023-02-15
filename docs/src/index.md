# HiQGA Documentation

This is a general purpose package for spatial statistical inference, geophysical forward modeling, Bayesian inference and inversion (both determinstic and probabilistic).

Readily usable geophysical forward operators are to do with AEM, CSEM and MT physics (references underneath), **for which the time domain AEM codes are fairly production-ready**. The current EM modeling is in 1D, but the inversion framework is dimensionally agnostic (e.g., you can regress images). Adding your own geophysical operators is easy, keep reading [down here](#Developing-HiQGA-or-modifying-it-for-your-own-special-forward-physics).

This package implements both the nested (2-layer) and vanilla trans-dimensional Gaussian process algorithm as published in 
- [*Bayesian inversion using nested trans-dimensional Gaussian processes*, A. Ray, Geophysical Journal International, **226(1)**, 2021](https://doi.org/10.1093/gji/ggab114).
- [*Bayesian geophysical inversion with trans-dimensional Gaussian process machine learning*, A. Ray and D. Myer, Geophysical Journal International **217(3)**, 2019](https://doi.org/10.1093/gji/ggz111).
- There is also a flavour of within-bounds Gauss-Newton/Occam's inversion implemented. For SkyTEM, TEMPEST and VTEM (all AEM), this is fully functional, but for other forward propagators you will have to provide a Jacobian (the linearization of the forward operator).

## Installation
[NCI](https://nci.org.au/) users look [here](#Development-setup-on-NCI) first!

To install, in a perfect world we'd use Julia's `Pkg` REPL by hitting `]` to enter `pkg>` mode. Then enter the following, at the `pkg>` prompt:
```
pkg> add HiQGA 
```
If installing to [follow along](https://www.youtube.com/watch?v=edgzr8vpCKY&list=PL0jP_ahe-BFmRWx6IT9G2zbFHA6qmJ52f&index=6) for the 2022 AEM [workshop branch](https://github.com/GeoscienceAustralia/HiQGA.jl/tree/workshop) use instead of the above `pkg> add HiQGA@0.2.2`

## Usage
Examples of how to use the package can be found in the `examples` directory. Simply `cd` to the relevant example directory and `include` the .`jl` files in the order they are named. If using VSCode make sure to do *Julia: Change to this Directory* from the three dots menu on the top right. The Markov Chain Monte Carlo sampler is configured to support parallel tempering on multiple CPUs - some of the examples accomplish this with Julia's built-in multiprocessing, and others use MPI in order to support inversions on HPC clusters that don't work with Julia's default SSH-based multiprocessing. The MPI examples require [MPI.jl](https://github.com/JuliaParallel/MPI.jl) and [MPIClusterManagers.jl](https://github.com/JuliaParallel/MPIClusterManagers.jl/), which are not installed as dependencies for this package, so you will need to ensure they are installed and configured correctly to run these examples. See [here](#Installing-MPI.jl-and-MPIClusterManagers.jl-on-NCI) for MPI on the NCI.

Some example scripts have as a dependency [Revise.jl](https://github.com/timholy/Revise.jl) as we're still actively [developing this package](https://pkgdocs.julialang.org/v1/getting-started/), so you may need to install [Revise](https://github.com/timholy/Revise.jl) if not already installed. All Julia users should be developing with [Revise](https://github.com/timholy/Revise.jl)! After installation, to run the examples, simply clone the package separately (or download as a ZIP), navigate to the `examples` folder and run the scripts in their numerical order.

## Updating the package 
```
pkg> update HiQGA
```

## Developing HiQGA or modifying it for your own special forward physics
After you have `]add`ed HiQGA you can simply do: 
```
pkg>dev HiQGA
```
If you haven't added it already, you can do:
```
pkg>dev https://github.com/GeoscienceAustralia/HiQGA.jl.git
```
It will download to your `JULIA_PKG_DEVDIR`. 

Another way is to simply clone or download this repository to your `JULIA_PKG_DEVDIR`, rename the cloned directory `HiQGA` removing the `.jl` bit and do
```
pkg>dev HiQGA
```
[Here's a gist](https://gist.github.com/a2ray/8c2c55c25fee6647501b403886bbe64d) on adding your own module if you want to modify the source code. Alternatively, if you only want to use the sampling methods in `HiQGA.transD_GP` without contributing to the source (boo! j/k) [here's another gist](https://gist.github.com/a2ray/92a8c14483c21dda6ddf56685b95fbb8) which is more appropriate. These gists were written originally for a package called `transD_GP` so you will have to modify `using transD_GP` to `using HiQGA.transD_GP`. Documentation is important and we're working on improving it before a full-release. 

## Development setup on NCI
You will need a Julia depot, where all packages are downloaded, and the package registry resides. While it may not be large in size, it can consume a lot of your inode (file count) quota. The easiest thing to do is set up a directory like this
```
mkdir/g/data/myprojectwithlotsofinodes/myusername/juliadepot
```
and then point a symlink to it from ***BOTH*** OOD and gadi, making sure you remove any existing `.julia` in your home directory with `rm -rf .julia` in your `$HOME`
```
cd
ln -s /g/data/myprojectwithlotsofinodes/myusername/juliadepot .julia
```
If you don't already have access to a `julia` binary, download the appropriate version `.tar.gz` from [here](https://julialang.org/downloads/) and then untar it in a location you have write access to. Then, in your `$HOME/bin` directory on **_BOTH_** OOD and gadi make a symlink to the julia binary like so:
```
cd ~/bin
ln -s /g/data/somwehere/julia-x.x.x/bin/julia .
```
The preferred development and usage environment for HiQGA is [Visual Studio Code](https://code.visualstudio.com/), which provides interactive execution of Julia code through the [VSCode Julia extension](https://code.visualstudio.com/docs/languages/julia). To install VSCode on the National Computational Infrastructure (NCI), you need to extract the VSCode rpm package using the steps in [this gist](https://gist.github.com/a2ray/701347f703b72abb630d2521b43c5f22), to a location where your account has write access. You will NOT be using vscode on a gadi login node, but on OOD.

Get Julia language support from VSCode after launching the VSCode binary by going to File->Extensions by searching for Julia. If after installation it doesn't find the Julia binary, go to File->Extensions->Julia->Manage(the little gear icon) and manually type in `/home/yourusername/bin/julia` in the "Executable Path" field.

It is also useful to use Revise.jl to ensure changes to the package are immediately reflected in a running Julia REPL (this is the reason that Revise is a dependency on some example scripts as noted above). More information on a workflow to use Revise during development can be found [here](https://gist.github.com/a2ray/e593751b24e45f8160ba8041fb811680).

**In your MPI job, make sure that you include in your qsub script** the `gdata` directory in which you have your julia executable and depot, e.g.,
```
#PBS -l storage=gdata/z67+gdata/kb5
```
### Installing MPI.jl and MPIClusterManagers.jl on NCI
We have found that the safest bet for MPI.jl to work without [UCX issues](https://docs.juliahub.com/MPI/nO0XF/0.19.2/knownissues/#UCX) on NCI is to use intel-mpi. In order to install MPI.jl and configure it to  use the intel-mpi provided by the module `intel-mpi/2019.8.254`, following the example below. 

```
$ module load intel-mpi/2019.8.254
$ julia

julia > ] 
pkg> add MPIPreferences
   Resolving package versions...
   Installed MPIPreferences ─ v0.1.7
    Updating `/g/data/up99/admin/yxs900/cr78_depot/environments/v1.7/Project.toml`
  [3da0fdf6] + MPIPreferences v0.1.7
    Updating `/g/data/up99/admin/yxs900/cr78_depot/environments/v1.7/Manifest.toml`
  [3da0fdf6] + MPIPreferences v0.1.7
Precompiling project...
  1 dependency successfully precompiled in 2 seconds (213 already precompiled)

julia > using MPIPreferences

julia> MPIPreferences.use_system_binary(;library_names=["/apps/intel-mpi/2019.8.254/intel64/lib/release/libmpi.so"],mpiexec="mpiexec",abi="MPICH",export_prefs=true,force=true)

┌ Info: MPI implementation identified
│   libmpi = "/apps/intel-mpi/2019.8.254/intel64/lib/release/libmpi.so"
│   version_string = "Intel(R) MPI Library 2019 Update 8 for Linux* OS\n"
│   impl = "IntelMPI"
│   version = v"2019.8.0"
└   abi = "MPICH"
┌ Info: MPIPreferences changed
│   binary = "system"
│   libmpi = "/apps/intel-mpi/2019.8.254/intel64/lib/release/libmpi.so"
│   abi = "MPICH"
└   mpiexec = "mpiexec"
```

Once the configuration is completed, install MPI.jl and MPIClusterManagers.jl.
```
pkg>add MPI, MPIClusterManagers, Distributed
```
Just to be safe, ensure that MPI has indeed built wth the version you have specified above:
```
julia> using MPI
julia> MPI.MPI_VERSION
v"3.1.0"

julia> MPI.MPI_LIBRARY
"IntelMPI"

julia> MPI.MPI_LIBRARY_VERSION
v"2019.8.0"

julia> MPI.identify_implementation()
("IntelMPI", v"2019.8.0")

```
To test, use an interactive NCI job with the following submission:
```
qsub -I -lwalltime=1:00:00,mem=16GB,ncpus=4,storage=gdata/z67+gdata/cr78
.
.
.
job is ready
```
now create a file called `mpitest.jl` with the following lines on some mount you have access to:
```
## MPI Init
using MPIClusterManagers, Distributed
import MPI
MPI.Init()
rank = MPI.Comm_rank(MPI.COMM_WORLD)
sz = MPI.Comm_size(MPI.COMM_WORLD)
if rank == 0
    @info "size is $sz"
end
manager = MPIClusterManagers.start_main_loop(MPI_TRANSPORT_ALL)
@info "there are $(nworkers()) workers"
MPIClusterManagers.stop_main_loop(manager)
rmprocs(workers())
exit()
```
Run the code after loading the intel-mpi module you have linked MPI.jl against with 
```
module load intel-mpi/2019.8.254
mpirun -np 3 julia mpitest.jl
```
and you should see output like:
```
[ Info: size is 3
[ Info: there are 2 workers
```
This is the basic recipe for all the cluster HiQGA jobs on NCI. After the call to `manager = MPIClusterManagers.start_main_loop(MPI_TRANSPORT_ALL)`, standard MPI execution stops, and we return to an explicit manager-worker mode with code execution only continuing on the manager which is Julia process 1.
### Installing PyPlot on NCI
Due to indode restrictions on NCI, we've resorted to using a communal matplotlib install as follows:
- Remove Conda, PyPlot, PyCall, HiQGA from your julia environment if it already exists
```
pkg> rm Conda
pkg> rm PyCall
pkg> rm PyPlot
pkg> rm HiQGA
```
- Delete the conda directory from your .julia directory (or wherever your julia depot is):
```
rm -rf conda/
```
- load python 3.8 on NCI and activate @richardt94 's virtual environment, then point julia at the python executable in this virtual env:
```
module load python3/3.8.5
source /g/data/z67/matplotlib-venv/bin/activate
PYTHON=/g/data/z67/matplotlib-venv/bin/python julia
```
Install and build PyCall:
```
pkg> add PyCall
pkg> build PyCall
julia> exit()
```
exit Julia and then restart Julia and in Pkg mode:
```
pkg> add PyPlot
```
- Install HiQGA in development mode:
```
pkg> dev HiQGA
```
### References for AEM and CSEM physics 

- [Blatter, D., Key, K., Ray, A., Foley, N., Tulaczyk, S., & Auken, E. (2018). Trans-dimensional Bayesian inversion of airborne transient EM data from Taylor Glacier, Antarctica. Geophysical Journal International, 214(3)](https://doi.org/10.1093/gji/ggy255)

- [Ray, A., & Key, K. (2012). Bayesian inversion of marine CSEM data with a trans-dimensional self parametrizing algorithm. Geophysical Journal International, 191(3), 1135-1151.](https://doi.org/10.1111/j.1365-246X.2012.05677.x)

