```@meta
CurrentModule = Microbiome
DocTestSetup  = quote
    using Microbiome
    using Microbiome.Dictionaries
end
```
# Samples and Features

Microbial profiles are made up of `AbstractSample`s and `AbstractFeature`s.
Typically, an `AbstractSample` is an individual biospecimen or other observation,
and contains some number of `AbstractFeature`s, such as taxa or gene functions.
`AbstractSample`s may also contain arbitrary metadata.

## Sample types

At its most basic, an `AbstractSample` simply encodes a `name`
(which should be a unique identifier), and a place to hold metadata.
The concrete type [`MicrobiomeSample`](@ref) is implemented with these two fields,
the latter of which is a `Dictionary` from [`Dictionaries.jl`](https://github.com/andyferris/Dictionaries.jl).

You can instantiate a `MicrobiomeSample` with just a name (in which case the metadata dictionary will be empty),
with existing metadata in the form of a dictionary,
or using keyword arguments for metadata entries.

```jldoctest sampletypes
julia> s1 = MicrobiomeSample("sample1")
MicrobiomeSample("sample1", {})

julia> s2 = MicrobiomeSample("sample2"; age=37)
MicrobiomeSample("sample2", {:age = 37})

julia> s3 = MicrobiomeSample("sample3", Dictionary([:gender, :age], ["female", 23]))
MicrobiomeSample("sample3", {:gender = "female", :age = 23})
```

### Working with metadata

To change or add metadata, you can use the [same syntax](https://github.com/andyferris/Dictionaries.jl#accessing-dictionaries)
as working with a [`Dictionary`] directly,
though note that this is a bit different from the `Dict` type in base julia:

```jldoctest sampletypes
julia> insert!(s1, :age, 50)
MicrobiomeSample("sample1", {:age = 50})

julia> set!(s3, :gender, "nonbinary")
MicrobiomeSample("sample3", {:gender = "nonbinary", :age = 23})

julia> delete!(s3, :gender)
MicrobiomeSample("sample3", {:age = 23})
```

You can access values of the dictionary using either `getindex`
or `getfield` syntax, that is:

```jldoctest sampletypes
julia> s3[:age]
23

julia> s3.age
23
```

Bulk addiction of metadata is also possible, by passing a `Dictionary` or `NamedTuple`
to `set!` or `insert!` (the latter will fail if any of the incoming keys are already found):

```jldoctest sampletypes
julia> insert!(s3, (gender = "nonbinary", genotype="XX"))
MicrobiomeSample("sample3", {:age = 23, :gender = "nonbinary", :genotype = "XX"})

julia> set!(s3, (genotype="XY", ses=7))
MicrobiomeSample("sample3", {:age = 23, :gender = "nonbinary", :genotype = "XY", :ses = 7})
```


## Feature Types

`AbstractFeature` types also have a `name`, but other fields are optional.
`Microbiome.jl` defines two concrete `AbstractFeature` types, [`Taxon`](@ref) and [`GeneFunction`](@ref).


### Taxon

The [`Taxon`](@ref) type contains a name and (optionally) a rank (eg `:phylum`).

```jldoctest taxon
julia> ecoli = Taxon("Escherichia_coli", :species)
Taxon("Escherichia_coli", :species)

julia> uncl = Taxon("Unknown_bug")
Taxon("Unknown_bug", missing)
```

You can access the name and rank fields using [`name`](@ref) and [`taxrank`](@ref) respectively, and also check whether the instance
has a rank with [`hasrank`](@ref), which returns `true` or `false`.

```jldoctest taxon
julia> hasrank(ecoli)
true

julia> hasrank(uncl)
false

julia> taxrank(ecoli)
:species

julia> taxrank(uncl)
missing
```

For compatibility with other tools, converting a `Taxon` to a `String`
will return the name prepended with the first letter of the taxonomic rank
and 2 underscores.
You can convert back using [`taxon`](@ref) (note the lowercase 't'):

```jldoctest taxon
julia> String(uncl)
"u__Unknown_bug"

julia> String(ecoli)
"s__Escherichia_coli"

julia> String(ecoli) |> Taxon
Taxon("s__Escherichia_coli", missing)

julia> String(ecoli) |> taxon
Taxon("Escherichia_coli", :species)
```



### GeneFunction

The [`GeneFunction`](@ref) type contains a name and (optionally) a [`Taxon`](@ref).
In addiction to providing both a name and `Taxon`,
you can instantiate a `GeneFunction` with just a name (in which case the taxon will be `missing`),
or with the name of the taxon (in which case it will not have a `rank`).

```jldoctest genefunction
julia> gf1 = GeneFunction("gene1")
GeneFunction("gene1", missing)

julia> gf2 = GeneFunction("gene2", "Species_name")
GeneFunction("gene2", Taxon("Species_name", missing))

julia> gf3 = GeneFunction("gene2", Taxon("Species_name", :species))
GeneFunction("gene2", Taxon("Species_name", :species))
```

You can access or check for various fields using similar methods as for `Taxon`:

```jldoctest genefunction
julia> hastaxon(gf1)
false

julia> hastaxon(gf2)
true

julia> hasrank(gf2)
false

julia> hasrank(gf3)
true

julia> name(gf3)
"gene2"

julia> taxon(gf3)
Taxon("Species_name", :species)

julia> taxrank(gf3)
:species
```

For compatibility with other tools,
Converting a `GeneFunction` to a `String` if it has a `Taxon`
will include the taxon name separated by `|`.
Converting back can be done using [`genefunction`](@ref)
(note the lowercase g and f).

```jldoctest genefunction
julia> String(gf3)
"gene2|s__Species_name"

julia> genefunction(String(gf3))
GeneFunction("gene2", Taxon("Species_name", :species))
```



## Types and Methods

```@docs
MicrobiomeSample
metadata
```

```@docs
Taxon
name
hasrank
taxrank
taxon
```

```@docs
GeneFunction
genefunction
```