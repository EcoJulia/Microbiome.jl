# Diversity measures on communities and samples

## Alpha Diversity

```@docs
shannon
shannon!
ginisimpson
ginisimpson!
```

## Beta Diversity

Quite often, it's useful to boil stuff down to distances between samples.
`AbstractAbundanceTable`s take advantage of the `pairwise` function
from [`Distances.jl`](https://github.com/JuliaStats/Distances.jl)
to get a symetric distance matrix.

Right now, only Bray-Curtis, Jaccard, and Hellinger are implemented, 
but it would be straightforward to add any others.
Open an issue if you want them!

You can also get fit a principal coordinates analysis (PCoA) to your `AbstractAbundanceTable`
using the `fit(MDS, ...)` from `MultivariateStats.jl` under the hood.

```@docs
braycurtis
jaccard
hellinger
pcoa
```
