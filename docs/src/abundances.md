## Working with microbial abundances

Tables of abundances are based off `ComMatrix` types from SpatialEcology.jl,
where columns are samples and rows are features (eg species). Sample and feature
names are also stored, and there's a convenience function if you want to convert
a `DataFrame` to a `ComMatrix`, assuming the first column contains feature
names:

```@example 1
using Microbiome
using DataFrames

df = DataFrame(species=["E. coli", "B. fragilis", "L. casei"],
                      sample1=[1, 4, 5],
                      sample2=[3, 8, 0],
                      sample3=[0, 3, 4])

abund = abundancetable(df)
```

Forgive the clutter... ComMatricies name rows as species (which is true in this
case, but need not be), and columns are "sites" rather than samples. That will
be fixed eventually.

```@example 1
@show samplenames(abund)

@show featurenames(abund)

@show sampletotals(abund) # this is column sums

@show featuretotals(abund)
```

If you want relative abundance, you can do `relativeabundance(abund)` or
`relativeabundance!(abund)`:

```@example 1
relativeabundance!(abund)

@show sampletotals(abund)
 ```

You can also filter on the `n` most abundant features accross the dataset. This
function automatically generates an `n+1` row for `other` containing the
remaining features. Note - these doesn't modify in-place, so you've gotta
reassign if you want to update:

```@example 1
abund2 = filterabund(abund, 1)

@show featurenames(abund2)
 ```

## Plotting

Some convenience plotting types are available using [`RecipesBase`][1] and
[StatPlots][2]

[1]: https://github.com/juliaplots/recipesbase.jl
[2]: https://github.com/juliaplots/StatPlots.jl

```@example 1
using StatPlots

abund = abundancetable(
           rand(100, 10), ["sample_$x" for x in 1:10],
           ["feature_$x" for x in 1:100])
relativeabundance!(abund)

abundanceplot(abund)

savefig("abundanceplot.png"); nothing # hide
```

![](abundanceplot.png)
