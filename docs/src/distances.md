# Working with Distances / Dissimilarity

Quite often, it's useful to boil stuff down to distances between samples.
`AbundanceTable`s can be used with the `pairwise()` function
from [`Distances.jl`](https://github.com/JuliaStats/Distances.jl)
to get a symetric distance matrix.

```@example 2
using Distances
using Microbiome

abund = abundancetable([1  3  0;
                        4  8  3;
                        5  0  4]);

dm = pairwise(BrayCurtis(), abund, dims=2)
```
