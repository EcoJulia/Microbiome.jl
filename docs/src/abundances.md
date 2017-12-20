## Working with microbial abundances

The `AbundanceTable` type is treated like a 2D array where columns are samples
and rows are features (eg species). Sample and feature names are also stored,
and there's a convenience function if you want to convert a `DataFrame` to an
`AbundanceTable`, assuming the first column contains feature names:

```julia
julia> using Microbiome
julia> using DataFrames

julia> df = DataFrame(species=["E. coli", "B. fragilis", "L. casei"],
                      sample1=[1, 4, 5],
                      sample2=[3, 8, 0],
                      sample3=[0, 3, 4])
3×4 DataFrames.DataFrame
│ Row │ species     │ sample1 │ sample2 │ sample3 │
├─────┼─────────────┼─────────┼─────────┼─────────┤
│ 1   │ E. coli     │ 1       │ 3       │ 0       │
│ 2   │ B. fragilis │ 4       │ 8       │ 3       │
│ 3   │ L. casei    │ 5       │ 0       │ 4       │

julia> abund = AbundanceTable(df)
3×3 Microbiome.AbundanceTable{Int64}:
 1  3  0
 4  8  3
 5  0  4
```

If you want relative abundance, you can do `relativeabundance(abund)` or
`relativeabundance!(abund)`:

```julia
julia> abund = relativeabundance(abund)
3×3 Microbiome.AbundanceTable{Float64}:
 0.1  0.272727  0.0
 0.4  0.727273  0.428571
 0.5  0.0       0.571429
 ```

You can also filter on the `n` most abundant features accross the dataset. This
function automatically generates an `n+1` row for `other` containing the
remaining features. Note - these doesn't modify in-place, so you've gotta
reassign if you want to update:

```julia
julia> abund2 = filterabund(abund, 1)
2×3 Microbiome.AbundanceTable{Float64}:
 0.4  0.727273  0.428571
 0.6  0.272727  0.571429

julia> abund2.features
2-element Array{String,1}:
 "B. fragilis"
 "other"
 ```

## Plotting

Some convenience plotting types are available using [`RecipesBase`][1] and
[StatPlots][2]

[1]: https://github.com/juliaplots/recipesbase.jl
[2]: https://github.com/juliaplots/StatPlots.jl

```julia
using StatPlots # TODO: add actual example

abund = AbundanceTable(
    rand(100, 10), ["sample_$x" for x in 1:10],
    ["feature_$x" for x in 1:100])

abund = relativeabundance(abund)
plot(abund, title="Random abundance")
```
