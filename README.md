[![Build Status](https://travis-ci.org/BioJulia/Microbiome.jl.svg?branch=master)](https://travis-ci.org/BioJulia/Microbiome.jl)

# Microbiome-related methods for julia

This package provides (or will provide) methods for microbial community
analyses. For now, I'm adding stuff as I need it, but pull-requests are more
than welcome.


## Working with Distances / Dissimilarity

Quite often, it's useful to boil stuff down to distances between samples. For
this, I'm using an interface with `Distances.jl` to generate a symetric
`DistanceMatrix`, which also contains a vector for samples, and a field
specifying which type of distance was used to calulate it. You can load one
in manually, or generate it from an `AbundanceTable`.

```julia
using Distances

dm = getdm(abund, BrayCurtis())
@show dm
@show dm.labels
```

I've also implemented a method to do a principle coordinates analysis. If
necessary, you can include `correct_neg=true` to use the correction method
described in [Lingoes (1971)][2]

```julia
p = pcoa(dm)

@show eigenvalue(p, 2)
@show principalcoord(p, 1)
```

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
