# HiQGA.jl

![CI status](https://github.com/GeoscienceAustralia/HiQGA.jl/workflows/CI/badge.svg)
[<img src="https://github.com/GeoscienceAustralia/HiQGA.jl/workflows/docs/badge.svg">](https://geoscienceaustralia.github.io/HiQGA.jl/)

This is a general purpose package for spatial statistical inference, geophysical forward modeling, Bayesian inference and inversion (both deterministic and probabilistic).

Readily usable geophysical forward operators are for AEM, CSEM and MT physics (references underneath), **for which the time domain AEM codes are fairly production-ready**. We've added [SMR](https://github.com/richardt94/SMRPInversion.jl) physics too! The current EM modeling is in 1D, but the inversion framework is dimensionally agnostic (e.g., you can regress images). Adding your own geophysical operators is [easy](https://geoscienceaustralia.github.io/HiQGA.jl/#Developing-HiQGA-or-modifying-it-for-your-own-special-forward-physics)! 

If you don't want to modify HiQGA at all, but simply add your own forward code for inference with HiQGA, check out the Shear Wave Dispersion [SWD](https://github.com/GeoscienceAustralia/HiQGA.jl/tree/master/examples/SWD) example that uses the `python` based [disba](https://pypi.org/project/disba/) library.

## Installation
To install the latest stable release, in a perfect world we'd use Julia's `Pkg` REPL by hitting `]` to enter `pkg>` mode. Then enter the following, at the `pkg>` prompt:
```
pkg> add HiQGA
```
Please note that as of v0.5.0, HiQGA will error for AEM inversions if the `electronics_halt.jl` is included with `read_survey_files()`. This is due to Julia language changes which disallow including from within a function (as the compiler has to deal with unexpected code changes during run-time -- probably a good thing). The error message tells you what to do if you used the old style `read_survey_files()`.

If you run into HiQGA build failures usually due to `PyPlot.jl` errors (`matplotlib` installation errors to do with `PyCall.jl` or `Conda.jl`), install `matplotlib` seperately with `Conda.jl`, before installing HiQGA and then you should be fine.

For the latest development version on here, you'd want to then do
```
pkg> dev HiQGA
```
## Docs
References, detailed instructions for installation, running examples and setting your environment on a cluster are ☞ [<img src="https://img.shields.io/badge/docs-stable-steelblue.svg">](https://geoscienceaustralia.github.io/HiQGA.jl/)

## Authors
- Anandaroop Ray
- Richard Taylor 

## Example AEM inversion
![](./aem.png)
