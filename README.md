# HiQGA.jl


**This branch is only for the 2022 Airborne Electromagnetics (AEM) workshop held in Perth as part of the [Exploring For the Future](https://www.eftf.ga.gov.au/) program [showcase](https://www.eftf.ga.gov.au/news/2022-showcase). See [this video](https://www.youtube.com/watch?v=edgzr8vpCKY&list=PL0jP_ahe-BFmRWx6IT9G2zbFHA6qmJ52f&index=6) using [these files](https://github.com/GeoscienceAustralia/HiQGA.jl/archive/refs/heads/workshop.zip) to run the examples.**

As part of this workshop, we'd recommend going through the following videos:
- An introduction to Geoscience Australia's (GA) ambitious 20 km spaced continent-wide AEM program by [Karol Czarnota](https://www.youtube.com/watch?v=pzJJf8RIipA&list=PL0jP_ahe-BFmRWx6IT9G2zbFHA6qmJ52f&index=1)
- How the Western Australian government has found 20 km spaced AEM data very useful, by [Klauss Gessner](https://www.youtube.com/watch?v=27YXK6RDkT0&list=PL0jP_ahe-BFmRWx6IT9G2zbFHA6qmJ52f&index=2)
- An introduction to AEM given by [Yusen Ley-Cooper](https://www.youtube.com/watch?v=KJxowEmCvHM&list=PL0jP_ahe-BFmRWx6IT9G2zbFHA6qmJ52f&index=3)
- An introduction to inverse theory presented by [Anandaroop Ray](https://www.youtube.com/watch?v=P2NhmWPQICQ&list=PL0jP_ahe-BFmRWx6IT9G2zbFHA6qmJ52f&index=5)
- Running the workshop Julia [examples](https://github.com/GeoscienceAustralia/HiQGA.jl/archive/refs/heads/workshop.zip) given in the `workshop` branch in the `examples` directory. See [here](https://www.youtube.com/watch?v=edgzr8vpCKY&list=PL0jP_ahe-BFmRWx6IT9G2zbFHA6qmJ52f&index=6) to follow along
- Integrate geophysics and geology in your subsurface interpretation as presented by [Sebastian Wong](https://www.youtube.com/watch?v=nsZ8IetMyew&list=PL0jP_ahe-BFmRWx6IT9G2zbFHA6qmJ52f&index=7)
- Avoid the 10 most common pitfalls in AEM interpretation according to [Neil Symington](https://www.youtube.com/watch?v=Of_-p6NIkJM&list=PL0jP_ahe-BFmRWx6IT9G2zbFHA6qmJ52f&index=8) 


## Installation
**Please use HiQGA v0.2.2 for this workshop!** Use Julia's `Pkg` REPL by hitting `]` to enter `pkg>` mode. Then enter the following, at the `pkg>` prompt:
```
pkg> add HiQGA@0.2.2
```

Unzip and run the Julia examples in the `examples` directory from the [zip file](https://github.com/GeoscienceAustralia/HiQGA.jl/archive/refs/heads/workshop.zip) of the `workshop` branch, by following the hands-on inversion video [here](https://www.youtube.com/watch?v=edgzr8vpCKY&list=PL0jP_ahe-BFmRWx6IT9G2zbFHA6qmJ52f&index=6) 

## Docs
References, detailed instructions for installation, running examples and setting your environment on a cluster are ☞ [<img src="https://img.shields.io/badge/docs-stable-steelblue.svg">](https://geoscienceaustralia.github.io/HiQGA.jl/)

## Example AEM inversion
![](./aem.png)
