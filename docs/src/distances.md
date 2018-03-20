# Working with Distances / Dissimilarity

Quite often, it's useful to boil stuff down to distances between samples. For
this, I'm using an interface with `Distances.jl` to generate a symetric
`DistanceMatrix`, which also contains a vector for samples, and a field
specifying which type of distance was used to calulate it. You can load one
in manually, or generate it from an `AbundanceTable`.

```julia
julia> using Distances
julia> using Microbiome

julia> abund = AbundanceTable([1  3  0;
                               4  8  3;
                               5  0  4])
3×3 Microbiome.AbundanceTable{Int64}:
 1  3  0
 4  8  3
 5  0  4

julia> dm = getdm(abund, BrayCurtis())
3×3 Microbiome.DistanceMatrix{Float64}:
 0.0       0.52381   0.176471
 0.52381   0.0       0.666667
 0.176471  0.666667  0.0
```

I've also implemented a method to do a principle coordinates analysis. If
necessary, you can include `correct_neg=true` to use the correction method
described in [Lingoes (1971)][2]

```julia
julia> p = pcoa(dm)

3×2 Microbiome.PCoA{Float64}:
 -0.251198   0.776895
  0.79841   -0.170903
 -0.547212  -0.605992

julia> eigenvalue(p, 2)
0.0050620487880063784

julia> principalcoord(p, 1)
3-element Array{Float64,1}:
 -0.251198
  0.79841
 -0.547212

julia> p.variance_explained
2-element Array{Float64,1}:
 0.979751
 0.0202492
```

## Plotting

Some convenience plotting types are available using [`RecipesBase`][1].

[1]: https://github.com/juliaplots/recipesbase.jl

```julia
using StatPlots # TODO: add actual example

abund = AbundanceTable(
    rand(100, 10), ["sample_$x" for x in 1:10],
    ["feature_$x" for x in 1:100])

dm = getdm(abund, BrayCurtis())
p = pcoa(dm, correct_neg=true)

plot(p, title="Random PCoA")
```

### Optimal Leaf Ordering

I've also provided a plotting recipe for making treeplots for `Hclust` objects
from the [`Clustering.jl`][2] package:

[2]: http://github.com/JuliaStats/Clustering.jl

```julia
dm = [0. .1 .2
      .1 0. .15
      .2 .15 0.]

h = hclust(dm, :single)
h.labels = ["a", "b", "c"]

plot(h)
```

![hclust plot 1](docs/img/hclustplot1.png)

Note that even though this is a valid tree, the leaf `a` is closer to leaf `c`,
despite the fact that `c` is more similar to `b` than to `a`. This can be fixed
with a method derived from the paper:

[Bar-Joseph et. al. "Fast optimal leaf ordering for hierarchical clustering." _Bioinformatics_. (2001)][3]

[3]: https://doi.org/10.1093/bioinformatics/17.suppl_1.S22

```julia
optimalorder!(h, dm)
plot(h)
```

![hclust plot 2](docs/img/hclustplot2.png)
