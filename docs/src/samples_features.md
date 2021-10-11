```@meta
CurrentModule=Microbiome
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

```jldoctest sample-types
julia> s1 = MicrobiomeSample("sample1")
MicrobiomeSample("sample1", {})

julia> s2 = MicrobiomeSample("sample2"; age=37)
MicrobiomeSample("sample2", {:age = 37})

julia> s3 = MicrobiomeSample("sample3", Dictionary([:gender, :age], ["female", 23]))
MicrobiomeSample("sample3", {:gender = "female", :age = 23})
```

To change or add metadata, you can use the [same syntax](https://github.com/andyferris/Dictionaries.jl#accessing-dictionaries)
as working with a [`Dictionary`] directly,
though note that this is a bit different from the `Dict` type in base julia:

```jldoctest sample-types
julia> insert!(s1, :age, 50)
MicrobiomeSample("sample1", {:age = 50})

julia> set!(s3, :gender, "nonbinary")
MicrobiomeSample("sample3", {:gender = "nonbinary", :age = 23})

julia> delete!(s3, :gender)
MicrobiomeSample("sample3", {:age = 23})
```


```@docs
MicrobiomeSample
```

## Feature Types

`AbstractFeature` types also have a `name`, but other fields are optional.
`Microbiome.jl` defines two concrete `AbstractFeature` types, [`Taxon`](@ref) and [`GeneFunction`](@ref).


### Taxon

The `Taxon` type contains a name and (optionally) a rank (eg `:phylum`).

```jldoctest
julia> ecoli = Taxon("Escherichia_coli", :species)
Taxon("Escherichia_coli", :species)

julia> uncl = Taxon("Unknown_bug")
Taxon("Unknown_bug", missing)
```

For compatibility with other tools, converting a `Taxon` to a `String`
will return the name prepended with the first letter of the taxonomic rank
and 2 underscores.
You can convert back using [`taxon`](@ref) (note the lowercase 't'):

```jldoctest
julia> String(uncl)
"u__Unknown_bug"

julia> String(ecoli)
"s__Escherichia_coli"

julia> String(ecoli) |> Taxon
Taxon("s__Escherichia_coli", missing)

julia> String(ecoli) |> taxon
Taxon("Escherichia_coli", :species)
```

```@docs
Taxon
GeneFunction
```

## Types and Methods

```@docs
metadata
name
rank
hasrank
taxon
hastaxon
```
