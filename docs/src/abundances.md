## Working with microbial abundances

Tables of abundances are based off `ComMatrix` types from SpatialEcology.jl,
where columns are samples and rows are features (eg species). Sample and feature
names are also stored, and there's a convenience function if you want to convert
a `DataFrame` to a `ComMatrix`, assuming the first column contains feature
names:

```@repl 1
using Microbiome
using DataFrames

df = DataFrame(species=["E. coli", "B. fragilis", "L. casei"],
                      sample1=[1, 4, 5],
                      sample2=[3, 8, 0],
                      sample3=[0, 3, 4]);

abund = abundancetable(df)
```

Forgive the clutter... ComMatricies name rows as species (which is true in this
case, but need not be), and columns are "sites" rather than samples. That will
be fixed eventually.

```@repl 1
samplenames(abund)
featurenames(abund)
sampletotals(abund) # this is column sums
featuretotals(abund) # this is row sums
```

If you want relative abundance, you can do `relativeabundance(abund)` or
`relativeabundance!(abund)`:

```@repl 1
relativeabundance!(abund)

sampletotals(abund)
```

You can also filter on the `n` most abundant features accross the dataset. This
function automatically generates an `n+1` row for `other` containing the
remaining features. Note - these doesn't modify in-place, so you've gotta
reassign if you want to update:

```@repl 1
abund2 = filterabund(abund, 1)

featurenames(abund2)
```

## Plotting

**NOTE: The following functions are not currently working - I've moved them to a new package to simplify dependencies. I'm leaving the docs for now as a reference - see `Microbiome.jl` versions 0.2.1 and below for working versions**

Some convenience plotting types are available using [`RecipesBase`](https://github.com/juliaplots/recipesbase.jl) and
[StatPlots](https://github.com/juliaplots/StatPlots.jl)

```@repl 1
using StatPlots

srand(1) # hide

abund = abundancetable(
           rand(100, 10),
           ["sample_$x" for x in 1:10],
           ["feature_$x" for x in 1:100]);
relativeabundance!(abund)

abundanceplot(abund)

savefig("abundanceplot.png"); nothing # hide
```

![](abundanceplot.png)
