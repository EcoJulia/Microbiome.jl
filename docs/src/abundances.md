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
                      sample3=[0, 3, 4]);

abund = abundancetable(df)
```

Forgive the clutter... ComMatricies name rows as species (which is true in this
case, but need not be), and columns are "sites" rather than samples. That will
be fixed eventually.

```@example 1
samplenames(abund)
featurenames(abund)
sampletotals(abund) # this is column sums
featuretotals(abund) # this is row sums
```

If you want relative abundance, you can do `relativeabundance(abund)` or
`relativeabundance!(abund)`:

```@example 1
relativeabundance!(abund)

sampletotals(abund)
```

You can also filter on the `n` most abundant features accross the dataset. This
function automatically generates an `n+1` row for `other` containing the
remaining features. Note - these doesn't modify in-place, so you've gotta
reassign if you want to update:

```@example 1
abund2 = filterabund(abund, 1)

featurenames(abund2)
```

## Plotting

Some convenience plotting types are available using
[MicrobiomePlots](https://github.com/BioJulia/MicrobiomePlots.jl)
and [StatsPlots](https://github.com/juliaplots/StatsPlots.jl)

```@example 1
ENV["GKSwstype"] = "100" # hide
using StatsPlots
using MicrobiomePlots
using Distributions
using Random # hide
Random.seed!(1) # hide

# add some high abundance bugs to be a bit more realistic
function spikein(spikes, y, x)
    m = rand(LogNormal(), y, x)
    for s in spikes
        m[s, :] = rand(LogNormal(3., 0.5), x)'
    end
    return m
end

# 100 species in 10 samples, with every 10th bug a bit more abundant
bugs = spikein(1:10:100, 100, 10)

abund = abundancetable(bugs,
    ["sample_$x" for x in 1:10],
    ["species_$x" for x in 1:100]);

relativeabundance!(abund)
abundanceplot(abund, xticks=(1:10, samplenames(abund)), xrotation=45)

savefig("abundanceplot.png"); nothing # hide
```

![abundance plot](./abundanceplot.png)

Perhaps you have some metadata that you'd like to add as well:

```@example 1
labels = ["a","a","b","a","b","b","b","b","a","a"]

plot(
    abundanceplot(abund, xticks=(1:10, samplenames(abund)), xrotation=45),
    plot(annotationbar(labels)),
    layout=grid(2,1, heights=[0.9,0.1]))

savefig("abundanceplot-annotations.png"); nothing # hide
```

![abundance plot with annotations](./abundanceplot-annotations.png)
