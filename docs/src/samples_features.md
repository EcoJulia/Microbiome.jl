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

```@docs
MicrobiomeSample
```

## Feature Types

`AbstractFeature` types also have a `name`, but other fields are optional.
`Microbiome.jl` defines two concrete `AbstractFeature` types, [`Taxon`](@ref) and [`GeneFunction`](@ref).

```@docs
Taxon
GeneFunction
```

## Types and Methods

```@docs
metadata
name
clade
hasclade
taxon
hastaxon
```
