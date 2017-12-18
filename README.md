[![Build Status](https://travis-ci.org/BioJulia/Microbiome.jl.svg?branch=master)](https://travis-ci.org/BioJulia/Microbiome.jl)

# Microbiome-related methods for julia

This package provides (or will provide) methods for microbial community
analyses. For now, I'm adding stuff as I need it, but pull-requests are more
than welcome.

## Plotting

I've included some plotting recipes for convenience using [`RecipesBase`][3].

```julia
using StatPlots

abund = AbundanceTable(
    rand(100, 10), ["sample_$x" for x in 1:10],
    ["feature_$x" for x in 1:100])

abund = relativeabundance(abund)
plot(abund, title="Random abundance")

dm = getdm(abund, BrayCurtis())
p = pcoa(dm, correct_neg=true)

plot(p, title="Random PCoA")
```



[1]: https://github.com/JuliaStats/Distances.jl/pull/76
[2]: https://doi.org/10.1007/BF02291398
[3]: https://github.com/juliaplots/recipesbase.jl
