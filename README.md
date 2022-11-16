# HiQGA.jl


**This branch of HiQGA is specifically for the 2022 Airborne Electromagnetics (AEM) workshop held in Perth as part of the [Exploring For the Future](https://www.eftf.ga.gov.au/) program [showcase](https://www.eftf.ga.gov.au/news/2022-showcase) held by [Geoscience Australia](https://www.ga.gov.au) (GA) in collaboration with the [Australian Institution of Geoscientists](https://www.aig.org.au/). See [this video](https://youtu.be/edgzr8vpCKY) using [these files](https://github.com/GeoscienceAustralia/HiQGA.jl/archive/refs/heads/workshop.zip) to run the examples.**

As part of this workshop, we'd recommend going through the following videos:
- Acknowledgement of country and an introduction to GA's 20 km spaced **continent-wide** AusAEM program by [Karol Czarnota](https://youtu.be/pzJJf8RIipA)
- How the Western Australia government has fruitfully used 20 km spaced AEM data, by [Klaus Gessner](https://youtu.be/27YXK6RDkT0)
- An introduction to AEM, surveying, and quality control given by [Yusen Ley-Cooper](https://youtu.be/KJxowEmCvHM)
- An introduction to inverse theory presented by [Anandaroop Ray](https://youtu.be/P2NhmWPQICQ)
- Running the workshop Julia [examples](https://github.com/GeoscienceAustralia/HiQGA.jl/archive/refs/heads/workshop.zip) given in the `workshop` branch in the `examples` directory. See [this video](https://youtu.be/edgzr8vpCKY) to follow along
- Integrate geophysics and geology in your subsurface interpretation as presented by [Sebastian Wong](https://youtu.be/nsZ8IetMyew)
- Avoid the 10 most common pitfalls in AEM interpretation according to [Neil Symington](https://youtu.be/Of_-p6NIkJM) 


## Installation
**Please use HiQGA v0.2.2 for this workshop!** Use Julia's `Pkg` REPL by hitting `]` to enter `pkg>` mode. Then enter the following, at the `pkg>` prompt:
```
pkg> add HiQGA@0.2.2
```

Unzip and run the Julia examples in the `examples` directory from the [zip file](https://github.com/GeoscienceAustralia/HiQGA.jl/archive/refs/heads/workshop.zip) of the `workshop` branch.

## Docs
References, detailed instructions for installation, running examples and setting your environment on a cluster are ☞ [<img src="https://img.shields.io/badge/docs-stable-steelblue.svg">](https://geoscienceaustralia.github.io/HiQGA.jl/). This workshop is indexed at GA as [eCat entry 147427](http://pid.geoscience.gov.au/dataset/ga/147427) and DOI [https://dx.doi.org/10.26186/147427](https://dx.doi.org/10.26186/147427)

## Example AEM inversion
![](./aem.png)
