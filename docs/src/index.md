# HiQGA Documentation

This is a general purpose package for spatial statistical inference, geophysical forward modeling, Bayesian inference and inversion (both determinstic and probabilistic).

Readily usable geophysical forward operators are to do with AEM, CSEM and MT physics (references underneath), **for which the time domain AEM codes are fairly production-ready**. The current EM modeling is in 1D, but the inversion framework is dimensionally agnostic (e.g., you e regress images). Adding your own geophysical operators is easy, keep reading [down here](#Developing-HiQGA-or-modifying-it-for-your-own-special-forward-physics).

This package implements both the nested (2-layer) and vanilla trans-dimensional Gaussian process algorithm as published in 
- [*An information theoretic Bayesian uncertainty analysis of AEM systems over Menindee Lakes, Australia*. A.Ray, A-Y. Ley-Cooper, R.C. Brodie, R. Taylor, N. Symington, N. F. Moghaddam, Geophysical Journal International, 235(2)., 2023](https://doi.org/10.1093/gji/ggad337).
- [*Bayesian inversion using nested trans-dimensional Gaussian processes*, A. Ray, Geophysical Journal International, **226(1)**, 2021](https://doi.org/10.1093/gji/ggab114).
- [*Bayesian geophysical inversion with trans-dimensional Gaussian process machine learning*, A. Ray and D. Myer, Geophysical Journal International **217(3)**, 2019](https://doi.org/10.1093/gji/ggz111).
- There is also a flavour of within-bounds Gauss-Newton/Occam's inversion implemented. For SkyTEM, TEMPEST and VTEM (all AEM), this is fully functional, but for other forward propagators you will have to provide a Jacobian (the linearization of the forward operator).

## Installation
Download Julia from [here](https://julialang.org/downloads/) and install the binary. HPC users e.g., on the [NCI](https://nci.org.au/) look [here, below in the document](#hpc-setup-on-nci) to set up MPI on a cluster for large jobs, or if you do not have Julia installed on your cluster.

Once you have started julia, to install HiQGA, use Julia's `Pkg` REPL by hitting `]` to enter `pkg>` mode within Julia like so: 
```
julia>]
(@v1.10) pkg>
```
Then enter the following, at the `pkg>` prompt:
```
pkg> add HiQGA 
```
If installing to [follow along](https://www.youtube.com/watch?v=edgzr8vpCKY&list=PL0jP_ahe-BFmRWx6IT9G2zbFHA6qmJ52f&index=6) for the 2022 AEM [workshop branch](https://github.com/GeoscienceAustralia/HiQGA.jl/tree/workshop) use instead of the above `pkg> add HiQGA@0.2.2`

## Usage
Examples of how to use the package can be found in the `examples` directory. Simply `cd` to the relevant example directory and `include` the .`jl` files in the order they are named. If using VSCode make sure to do *Julia: Change to this Directory* from the three dots menu on the top right. The Markov Chain Monte Carlo sampler is configured to support parallel tempering on multiple CPUs - some of the examples accomplish this with Julia's built-in multiprocessing, and others use MPI in order to support inversions on HPC clusters that don't work with Julia's default SSH-based multiprocessing. The MPI examples require [MPI.jl](https://github.com/JuliaParallel/MPI.jl) and [MPIClusterManagers.jl](https://github.com/JuliaParallel/MPIClusterManagers.jl/), which are not installed as dependencies for this package, so you will need to ensure they are installed and configured correctly to run these examples. See [here, below in the document](#installing-mpijl-and-mpiclustermanagersjl-on-nci) for MPI on the NCI.

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

## HPC setup on NCI
If you don't already have access to a `julia` binary, download the appropriate version `.tar.gz` from [here](https://julialang.org/downloads/) and then untar it in a location you have write access to, like so: 
```
cd /somwehere/home/me
tar -xvzf /somwehere/home/me/julia-x.x.x.tar.gz
```
Then, in your `$HOME/bin` directory make a symlink to the julia binary like so:
```
cd ~/bin
ln -s /somwehere/home/me/julia-x.x.x/bin/julia .
```
Make sure your `$HOME/bin` is in your `$PATH` else which you can check with `echo $PATH | grep "$HOME/bin"`. If you do not see your `bin` directory highlighted, do `export PATH=~/bin:$PATH`
The preferred development and usage environment for HiQGA is [Visual Studio Code](https://code.visualstudio.com/), which provides interactive execution of Julia code through the [VSCode Julia extension](https://code.visualstudio.com/docs/languages/julia). To install VSCode on the National Computational Infrastructure (NCI), you need to extract the VSCode rpm package using the steps in [this gist](https://gist.github.com/a2ray/701347f703b72abb630d2521b43c5f22), to a location where your account has write access. You will NOT be using vscode on a gadi login node, but on OOD.

Get Julia language support from VSCode after launching the VSCode binary by going to File->Extensions by searching for Julia. If after installation it doesn't find the Julia binary, go to File->Extensions->Julia->Manage (the little gear icon) and manually type in `/home/yourusername/bin/julia` in the "Executable Path" field.

It is also useful to use Revise.jl to ensure changes to the package are immediately reflected in a running Julia REPL (this is the reason that Revise is a dependency on some example scripts as noted above). More information on a workflow to use Revise during development can be found [here](https://gist.github.com/a2ray/e593751b24e45f8160ba8041fb811680).

As [shown above](#installation), to install HiQGA, start Julia and go into the `Pkg` REPL by hitting `]` to enter `pkg>` mode like so: 
```
julia>]
(@v1.10) pkg>
```
Then enter the following, at the `pkg>` prompt:
```
pkg> add HiQGA 
```

### Installing MPI.jl and MPIClusterManagers.jl on NCI
We have found that the safest bet for MPI.jl to work without [UCX issues](https://docs.juliahub.com/MPI/nO0XF/0.19.2/knownissues/#UCX) on NCI is to use intel-mpi. In order to install MPI.jl and configure it to  use the intel-mpi provided by the module `intel-mpi/2021.10.0`, following the example below. 

```
$ module load intel-mpi/2021.10.0
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

julia> MPIPreferences.use_system_binary(;library_names=["/apps/intel-mpi/2021.10.0/lib/release/libmpi.so"],mpiexec="mpiexec",abi="MPICH",export_prefs=true,force=true)
[Info: MPIPreferences changed
|   binary = "system"
|   libmpi = "/apps/intel-mpi/2021.10.0/lib/release/libmpi.so"
|   abi = "MPICH"
|   mpiexec = "mpiexec"
|   preloads = Any[]
[   preloads_env_switch = nothing

julia> exit()
```
Once the configuration is completed, install MPI.jl and MPIClusterManagers.jl in a restarted Julia session. We had errors with other versions of MPI.jl besides v0.19.2 on NCI, maybe not an issue elsewhere.
```
pkg>add MPI@0.19.2, MPIClusterManagers, Distributed
Resolving package versions...
  No Changes to `~/.julia/environments/v1.9/Project.toml`
  No Changes to `~/.julia/environments/v1.9/Manifest.toml`
Precompiling project...
  4 dependencies successfully precompiled in 5 seconds. 225 already precompiled.
```
Just to be safe, ensure that MPI has indeed built wth the version you have specified above:
```
julia> using MPI

julia> MPI.MPI_VERSION
v"3.1.0"

julia> MPI.MPI_LIBRARY
IntelMPI::MPIImpl = 4

julia> MPI.MPI_LIBRARY_VERSION
v"2021.0.0"

julia> MPI.identify_implementation()
(MPI.IntelMPI, v"2021.0.0")

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
module load intel-mpi/2021.10.0
mpirun -np 3 julia mpitest.jl
```
and you should see output like:
```
[ Info: size is 3
[ Info: there are 2 workers
```
This is the basic recipe for all the cluster HiQGA jobs on NCI. After the call to `manager = MPIClusterManagers.start_main_loop(MPI_TRANSPORT_ALL)`, standard MPI execution stops, and we return to an explicit manager-worker mode with code execution only continuing on the manager which is Julia process 1.

The next time you start julia you have HiQGA ready for use with
```
julia> using HiQGA
```
navigate to the [examples](https://github.com/GeoscienceAustralia/HiQGA.jl/tree/master/examples) folder to run some example scripts. **You can end here as a regular user, however for development mode see below.**

### MPI or parallel considerations for number of CPUs
For the deterministic AEM inversions, you can do as many inversions in parallel as you have worker nodes. With McMC though, and the parallel tempering, it is a little more complicated. First of all, if you have a node with 48 CPUs (i.e., ppn=48), it is best to keep an integer number of soundings within a node. Helicopter systems typically need 5 chains per sounding (or, 1 manager + 5 workers = 6 CPUs/sounding) and fixed wing systems need 7 chains per sounding (or, 1 manager + 7 workers = 8 CPUs/sounding) as the likelihood for Tx-Rx geometry misfits can be a little rugose. Secondly, the variable `nchainspersounding+1` must exactly divide the total number of available CPU **workers** when using the function `loopacrossAEMsoundings`. From any numberof cpus total provided to HiQGA, processes `1` through `1+nchainspersounding` are reserved for memory purposes on for massive production jobs. Therefore, to see what works, try running the following before submitting or starting a parallel job:
```julia
using HiQGA
# say there are 24 CPUs on a large compute box you can use
# say you've done this with an addprocs(23), then for 3 soundings and a helicopter system,
transD_GP.splittasks(;nsoundings=3, ncores=23, nchainspersounding=5, ppn=24)
```
You will get
```julia
[ Info: will require 1 iterations of 3 soundings in one iteration
```
Say for an actual MPI job where you have a thousands of cores for a fixed wing system, say `ppn=104` and you have 10 nodes available with 1040 cpus. Run this before starting, to see how many iterations it will take to do 800 soundings.
```julia
using HiQGA
transD_GP.splittasks(;nsoundings=800, ncores=1039, nchainspersounding=7, ppn=104)
[ Info: will require 7 iterations of 129 soundings in one iteration
```
You can also violate the stay within a node principle if you have a cluster that does efficient message passing, like this. Say you have an actual `ppn=104` as in the case of the NCI Sandy Rapids nodes. Then for 12 nodes = 104*12 = 1248 CPUs, we ensure this total number is divisible by the `ppn` we are providing HiQGA. Say this is 48, a nice number for helicopter systems since 48/(5+1) is an integer. We can fool HiQGA into thinking we are staying within nodes of 48 CPUs as it only uses the `ppn` value we provide it. We can request 1248 CPUs through `PBS` and a `qsub` script but check our input first like so:
```julia
transD_GP.splittasks(;nsoundings=1211, ncores=1247, nchainspersounding=5, ppn=48)
[ Info: will require 6 iterations of 207 soundings in one iteration
```
[Here](https://github.com/GeoscienceAustralia/HiQGA.jl/blob/a8b258d6cef23be7423c9e8652ea0926af28f448/ASEG_Hobart_Workshop_2024/UDF_probabilistic/submit.sh) is an example of a massive job qsub submit script.
### Troubleshooting MPI set up
Some folks have reported that the above MPI install provides error messages to the order of "You are using the system provided MPI", indicating that it is not Intel MPI 2021.10.0 that they are working with. In this case, you should first remove MPI.jl
```
# goto Pkg mode in Julia by hitting ]
pgk> rm MPI.jl
```
exit Julia, then edit `~/.julia/prefs/MPI.toml` adding in the following lines
```
path = "/apps/intel-mpi/2021.10.0"
library = "/apps/intel-mpi/2021.10.0/lib/release/libmpi.so"
binary = "system"
```
then go back to Julia in Pkg mode, making sure the intel-mpi module is loaded in BASH and add back MPI.jl
```
module load intel-mpi/2021.10.0
julia
```
Within Julia, you can then do 
```
# hit ] to enter Pkg mode
pkg>add MPI@0.19.2
```
and this should work.
## For installing development mode pre-release versions
```
pkg> dev HiQGA
```
**Make a pull request if you wish to submit your change -- we welcome feature additions**. If you want to switch back to the official version from development mode, do
```
pkg> free HiQGA
```
## References for AEM and CSEM physics 

- [Blatter, D., Key, K., Ray, A., Foley, N., Tulaczyk, S., & Auken, E. (2018). Trans-dimensional Bayesian inversion of airborne transient EM data from Taylor Glacier, Antarctica. Geophysical Journal International, 214(3)](https://doi.org/10.1093/gji/ggy255)
- [Brodie, R. C. (2010). Holistic inversion of airborne electromagnetic data. PhD thesis, Australian National University](https://openresearch-repository.anu.edu.au/bitstream/1885/49403/4/02Whole.pdf).
- [Ray, A., & Key, K. (2012). Bayesian inversion of marine CSEM data with a trans-dimensional self parametrizing algorithm. Geophysical Journal International, 191(3), 1135-1151.](https://doi.org/10.1111/j.1365-246X.2012.05677.x)

