# Working with Distances / Dissimilarity

Quite often, it's useful to boil stuff down to distances between samples. For
this, I'm using an interface with `Distances.jl` to generate a symetric
`DistanceMatrix`, which also contains a vector for samples, and a field
specifying which type of distance was used to calulate it. You can load one
in manually, or generate it from an `AbundanceTable`.

```@repl 2
using Distances
using Microbiome

abund = abundancetable([1  3  0;
                        4  8  3;
                        5  0  4]);

dm = getdm(abund, BrayCurtis())
```

I've also implemented a method to do a principle coordinates analysis. If
necessary, you can include `correct_neg=true` to use the correction method
described in [Lingoes (1971)](http://dx.doi.org/10.1007/BF02291398)

```@repl 2
p = pcoa(dm)

eigenvalue(p, 2)
principalcoord(p, 1)
variance(p, [1,2])
```

## Plotting

**NOTE: The following functions are not currently working - I've moved them to a new package to simplify dependencies. I'm leaving the docs for now as a reference - see `Microbiome.jl` versions 0.2.1 and below for working versions**

Some convenience plotting types are available using [`RecipesBase`](https://github.com/juliaplots/recipesbase.jl).

```@repl 2
using StatPlots

srand(1) # hide
abund = abundancetable(
    rand(100, 10),
    ["sample_$x" for x in 1:10],
    ["feature_$x" for x in 1:100]);

dm = getdm(abund, BrayCurtis());
p = pcoa(dm, correct_neg=true);

plot(p, title="Random PCoA")
savefig("pcoplot.png"); nothing # hide
```

![pcoa plot](pcoplot.png)

### Optimal Leaf Ordering

I've also provided a plotting recipe for making treeplots for `Hclust` objects
from the [`Clustering.jl`](http://github.com/JuliaStats/Clustering.jl) package:

```@repl 2
using Clustering

dm = [0. .1 .2
      .1 0. .15
      .2 .15 0.];

h = hclust(dm, :single);
h.labels = ["a", "b", "c"];

hclustplot(h)
savefig("hclustplot1.png"); nothing # hide
```

![hclust plot 1](hclustplot1.png)

Note that even though this is a valid tree, the leaf `a` is closer to leaf `c`,
despite the fact that `c` is more similar to `b` than to `a`. This can be fixed
with a method derived from the paper:

[Bar-Joseph et. al. "Fast optimal leaf ordering for hierarchical clustering." _Bioinformatics_. (2001)](https://doi.org/10.1093/bioinformatics/17.suppl_1.S22)

```@repl 2
optimalorder!(h, dm)
hclustplot(h)

savefig("hclustplot2.png"); nothing # hide
```

![hclust plot 1](hclustplot2.png)
