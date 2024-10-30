# [HiQGA.jl](https://github.com/GeoscienceAustralia/HiQGA.jl) workshop at ASEG Discover 2024 conference, Hobart: Airborne electromagnetic (AEM) surveying from farm to table -- an interactive overview of acquisition, processing and inversion of AEM data

Jupyter notebooks are in the directory `00_notebooks/` and will only require the use of a computer with at least 7 CPUs. Most laptops can run these in 2024. Remember to change any file paths to your local disk if running the notebooks yourself. To view the notebooks, simply click the `.ipynb` file in GitHub on your browser.

Julia scripts on which the notebooks are based, are in the directories

```
synthetic/  
UDF_data/  
UDF_deterministic/  
UDF_probabilistic/  
```

Simply run the files in numbered order for `synthetic/` and `UDF_data/`. Remember to change any file paths to your local disk!

Please note that the `UDF_deterministic/` and `UDF_probabilistic/` directories use MPI and 1,000s of cores. The submit scripts for qsub/PBS on the Gadi cluster, as well as diagnostic output from HiQGA is included. These two directories only require the 01, 02 and 04 scripts to plot pre-run results. However, for the probabilistic inversions, the full posteriors are too large to provide here. Again, remember to change any file paths to your local disk.

## Installation of HiQGA version used in the workshop
To install the latest stable release, in a perfect world we'd use Julia's `Pkg` REPL by hitting `]` to enter `pkg>` mode. Then enter the following, at the `pkg>` prompt:
```
pkg> add HiQGA @0.4.8
```
For the latest HiQGA version (not used in the workshop), drop the `@0.4.8`

## Docs
References, detailed instructions for installation, running examples and setting your environment on a cluster are ☞ [<img src="https://img.shields.io/badge/docs-stable-steelblue.svg">](https://geoscienceaustralia.github.io/HiQGA.jl/)

## References
Please see the accompanying conference abstract for mathmatical details and further references
https://doi.org/10.5281/zenodo.13918176



