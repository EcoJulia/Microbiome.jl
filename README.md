[![Build Status](https://travis-ci.org/kescobo/Microbiome.jl.svg?branch=master)](https://travis-ci.org/kescobo/Microbiome.jl)

# Microbiome-related methods for julia

This package provides (or will provide) methods for microbial community
analyses. For now, I'm adding stuff as I need it, but pull-requests are more
than welcome.

## Working with microbial abundances

The `AbundanceTable` type is treated like a 2D array where columns are samples
and rows are features (eg species). Sample and feature names are also stored,
and there's a convenience function if you want to convert a `DataFrame` to an
`AbundanceTable`, assuming the first column contains feature names:

```julia
using Microbiome
using DataFrames

df = DataFrame(species=["E. coli", "B. fragilis", "L. casei"],
               sample1=[0.1, 0.4, 0.5],
               sample2=[0.3, 0.7, 0.0],
               sample3=[0.0, 0.2, 0.8])

abund = abundancetable(df)

@show abund
full(abund)
@show sitenames(abund)
@show specnames(features)
```

Note: I've used a relative abundance table here, but that need not be the case.
You can also use raw counts, RPK etc. If you want relative abundance, you can
do `relativeabundance(abund)`

You can also filter on the `n` most abundant features accross the dataset. This
function automatically generates an `n+1` row for `other` containing the
remaining features. Note - these doesn't modify in-place, so you've gotta
reassign if you want to update:

```julia
abund2 = filterabund(abund, 1)
@show abund2
full(abund2)
@show specnames(abund2)
```

## Working with Distances / Dissimilarity

Quite often, it's useful to boil stuff down to distances between samples. For
this, I'm using an interface with `Distances.jl` to generate a symetric
`DistanceMatrix`, which also contains a vector for samples, and a field
specifying which type of distance was used to calulate it. You can load one
in manually, or generate it from an `AbundanceTable` (note: `BrayCurtis` depends
on a [PR][1] that's not yet merged into `Distances.jl`)

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

abund = abundancetable(rand(100, 10))

abund = relativeabundance(abund)
plot(abund, title="Random abundance")

dm = getdm(abund, BrayCurtis())
p = pcoa(dm, correct_neg=true)

plot(p, title="Random PCoA")
```



[1]: https://github.com/JuliaStats/Distances.jl/pull/76
[2]: https://doi.org/10.1007/BF02291398
[3]: https://github.com/juliaplots/recipesbase.jl
