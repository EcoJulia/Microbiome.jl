```@meta
CurrentModule = Microbiome
DocTestSetup  = quote
    using Microbiome
    using Microbiome.SparseArrays
    using Random
    Random.seed!(42)
end
```

# Working with microbial abundances

The primary type for working with microbial abundances is the [`CommunityProfile`](@ref),
which is a sparse matrix with [`MicrobiomeSample`](@ref)s as column indices
and features (eg [`Taxon`](@ref)s or [`GeneFunction`](@ref)s) as row indices.

For example, let's make a `CommunityProfile` with 3 samples,
5 species and 5 genera.

```jldoctest profiles
julia> samps = MicrobiomeSample.(["s1", "s2", "s3"])
3-element Vector{MicrobiomeSample}:
 MicrobiomeSample("s1", {})
 MicrobiomeSample("s2", {})
 MicrobiomeSample("s3", {})

julia> taxa = [[Taxon("s$i", :species) for i in 1:5]; [Taxon("g$i", :genus) for i in 1:5]]
10-element Vector{Taxon}:
 Taxon("s1", :species)
 Taxon("s2", :species)
 Taxon("s3", :species)
 Taxon("s4", :species)
 Taxon("s5", :species)
 Taxon("g1", :genus)
 Taxon("g2", :genus)
 Taxon("g3", :genus)
 Taxon("g4", :genus)
 Taxon("g5", :genus)

julia> mat = spzeros(10, 3); # 10 x 3 matrix filled with zeros

julia> for i in 1:10, j in 1:3 
           # fill some spots with random values
           rand() < 0.3 && (mat[i,j] = rand())
       end

julia> mat
10×3 SparseMatrixCSC{Float64, Int64} with 10 stored entries:
  ⋅         ⋅        0.172933
  ⋅         ⋅         ⋅
 0.956916   ⋅         ⋅
 0.422956   ⋅         ⋅
  ⋅         ⋅         ⋅
 0.502952   ⋅        0.167169
  ⋅         ⋅        0.244683
  ⋅         ⋅        0.143638
 0.570085  0.249238   ⋅
  ⋅        0.841643   ⋅

julia> comm = CommunityProfile(mat, taxa, samps)
CommunityProfile{Float64, Taxon, MicrobiomeSample} with 10 features in 3 samples

Feature names:
s1, s2, s3...g4, g5

Sample names:
s1, s2, s3
```

## Accessing `CommunitProfile` contents

It is easy to get out the underlying abundance matrix,
features, and samples using
[`abundances`](@ref), [`features`](@ref), and [`samples`](@ref) respectively:

```jldoctest profiles
julia> abundances(comm)
10×3 SparseMatrixCSC{Float64, Int64} with 10 stored entries:
  ⋅         ⋅        0.172933
  ⋅         ⋅         ⋅
 0.956916   ⋅         ⋅
 0.422956   ⋅         ⋅
  ⋅         ⋅         ⋅
 0.502952   ⋅        0.167169
  ⋅         ⋅        0.244683
  ⋅         ⋅        0.143638
 0.570085  0.249238   ⋅
  ⋅        0.841643   ⋅

julia> features(comm)
10-element Vector{Taxon}:
 Taxon("s1", :species)
 Taxon("s2", :species)
 Taxon("s3", :species)
 Taxon("s4", :species)
 Taxon("s5", :species)
 Taxon("g1", :genus)
 Taxon("g2", :genus)
 Taxon("g3", :genus)
 Taxon("g4", :genus)
 Taxon("g5", :genus)
 
julia> samples(comm)
3-element Vector{MicrobiomeSample}:
 MicrobiomeSample("s1", {})
 MicrobiomeSample("s2", {})
 MicrobiomeSample("s3", {})
```

You can also just get an array of feature and sample names
using [`featurenames`](@ref) and [`samplenames`](@ref) respectively.

```jldoctest profiles
julia> featurenames(comm)
10-element Vector{String}:
 "s1"
 "s2"
 "s3"
 "s4"
 "s5"
 "g1"
 "g2"
 "g3"
 "g4"
 "g5"

julia> samplenames(comm)
3-element Vector{String}:
 "s1"
 "s2"
 "s3"
```

Finally, you can pull out the metadata of all samples
or a subset using [`metadata`](@ref).
The returned value is a vector of `NamedTuple`s,
which is compliant with the [`Tables.jl`](https://github.com/JuliaData/Tables.jl) interface,
so it's easy to load into other formats (like [`DataFrames.jl`](https://github.com/JuliaData/DataFrames.jl) for example):

```jldoctest profiles
julia> metadata(comm)
3-element Vector{NamedTuple{(:sample,), Tuple{String}}}:
 (sample = "s1",)
 (sample = "s2",)
 (sample = "s3",)
```

Of course, for now, the metadata is pretty boring - 
jump ahead to [Working with metadata](@ref working-metadata)
to do some more fun things.

## Indexing and selecting

`CommunityProfile`s wrap a sparse matrix,
and you can access the values as you would a normal matrix.
In julia, you can pull out specific values using `[row, col]`.
So for example, to get the 3rd row, 2nd column, of matrix `mat`:

```jldoctest indexing
julia> mat = reshape(1:12, 4, 3) |> collect
4×3 Matrix{Int64}:
 1  5   9
 2  6  10
 3  7  11
 4  8  12

julia> mat[3,2]
7
```

You can also get "slices", eg to get rows 2-4, column 1:

```jldoctest indexing
julia> mat[2:4, 1]
3-element Vector{Int64}:
 2
 3
 4
```

To get all of one dimension, you can just use a bare `:`

```jldoctest indexing
julia> mat[:, 1:2]
4×2 Matrix{Int64}:
 1  5
 2  6
 3  7
 4  8
```

For `CommunityProfile`s, indexing with integer values
will return the value of the matrix at that position,
while indexing with slices will return a new `CommunityProfile`:

```jldoctest profiles
julia> comm[1,3]
0.17293302893695128

julia> comm[6:8,3]
CommunityProfile{Float64, Taxon, MicrobiomeSample} with 3 features in 1 samples

Feature names:
g1, g2, g3

Sample names:
s3



julia> comm[1:3,3] |> abundances
3×1 SparseMatrixCSC{Float64, Int64} with 1 stored entry:
 0.172933
  ⋅
  ⋅
```

### Indexing with strings and regular expressions

It is often inconvenient to find the numerical index
of a particular feature or sample.
Instead, you can use strings or regular expressions
to get slices of a `CommunityProfile`,
which will match on the `name` field of the features or samples.
This kind of indexing always returns a `CommunityProfile`,
even if it only has 1 value.

```jldoctest profiles
julia> comm["g1", "s1"]
CommunityProfile{Float64, Taxon, MicrobiomeSample} with 1 features in 1 samples

Feature names:
g1

Sample names:
s1



julia> comm[r"[gs]1", "s1"]
CommunityProfile{Float64, Taxon, MicrobiomeSample} with 2 features in 1 samples

Feature names:
s1, g1

Sample names:
s1
```



## [Working with metadata](@id working-metadata)

The `metadata` dictionaries of `MicrobiomeSample`s can accessed
and updated inside a `CommunityProfile` using `insert!`, `delete!`, `set!`, and `unset!`.
either by directly acting on the underlying sample object,
or using special methods that take the sample name as the second argument:

```jldoctest profiles
julia> set!(samples(comm)[1], :subject, "kevin")
MicrobiomeSample("s1", {:subject = "kevin"})

julia> insert!(comm, "s2", :subject, "anika")
MicrobiomeSample("s2", {:subject = "anika"})
```

You can also retrieve all metadata associated with a table using [`metadata`](@ref),
which will return a `Tables.jl`-compatible vector or `NamedTuple`s,
where each "row" corresponds to one sample.
All metadata fields found in any sample will be returned in every row,
with the value `missing` in any samples that do not have that field set.

```jldoctest profiles
julia> metadata(comm)
3-element Vector{NamedTuple{(:sample, :subject), T} where T<:Tuple}:
 (sample = "s1", subject = "kevin")
 (sample = "s2", subject = "anika")
 (sample = "s3", subject = missing)
```

And you can bulk-`insert!` or `set!` metadata by passing a similar Table-like object
with the a field (`:sample` by default) matching sample names found in the `CommunityProfile`.

```jldoctest profiles
julia> md = [(sample="s1", subject="kevin", foo="bar"), (sample="s3", subject="annelle", foo="baz")]
2-element Vector{NamedTuple{(:sample, :subject, :foo), Tuple{String, String, String}}}:
 (sample = "s1", subject = "kevin", foo = "bar")
 (sample = "s3", subject = "annelle", foo = "baz")

julia> set!(comm, md)

julia> md2 = [(name="s1", other="Hello, World!"), (name="s2", other="Goodbye!")]
2-element Vector{NamedTuple{(:name, :other), Tuple{String, String}}}:
 (name = "s1", other = "Hello, World!")
 (name = "s2", other = "Goodbye!")

julia> insert!(comm, md2; namecol=:name)

julia> metadata(comm)
3-element Vector{NamedTuple{(:sample, :subject, :foo, :other), T} where T<:Tuple}:
 (sample = "s1", subject = "kevin", foo = "bar", other = "Hello, World!")
 (sample = "s2", subject = "anika", foo = missing, other = "Goodbye!")
 (sample = "s3", subject = "annelle", foo = "baz", other = missing)
```

## Types and Methods

```@docs
CommunityProfile
samples
features
samplenames
featurenames
commjoin
relativeabundance
relativeabundance!
present
prevalence
prevalence_filter
```
