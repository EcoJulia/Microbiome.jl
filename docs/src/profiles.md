```@meta
CurrentModule = Microbiome
DocTestSetup  = quote
    using Microbiome
    using Microbiome.SparseArrays

    open("taxprof.csv", "w") do io
        print(io, """
        features,s1,s2
        s__Escherichia_coli,1,0
        s__Bifidobacterium_longum,3,2
        s__Prevotella_copri,5,4
        """)
    end
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

julia> using SparseArrays

julia> mat = spzeros(10, 3); # 10 x 3 matrix filled with zeros

julia> for i in 1:10, j in 1:3
           if all(isodd, (i,j))
               mat[i,j] = i+j
           elseif all(iseven, (i, j))
               mat[i,j] = i / j
           end
       end

julia> mat
10×3 SparseMatrixCSC{Float64, Int64} with 15 stored entries:
  2.0   ⋅    4.0
   ⋅   1.0    ⋅
  4.0   ⋅    6.0
   ⋅   2.0    ⋅
  6.0   ⋅    8.0
   ⋅   3.0    ⋅
  8.0   ⋅   10.0
   ⋅   4.0    ⋅
 10.0   ⋅   12.0
   ⋅   5.0    ⋅

julia> comm = CommunityProfile(mat, taxa, samps)
CommunityProfile{Float64, Taxon, MicrobiomeSample} with 10 features in 3 samples

Feature names:
s1, s2, s3...g4, g5

Sample names:
s1, s2, s3
```

## Accessing `CommunityProfile` contents

It is easy to get out the underlying abundance matrix,
features, and samples using
[`abundances`](@ref), [`features`](@ref), and [`samples`](@ref) respectively:

```jldoctest profiles
julia> abundances(comm)
10×3 SparseMatrixCSC{Float64, Int64} with 15 stored entries:
  2.0   ⋅    4.0
   ⋅   1.0    ⋅
  4.0   ⋅    6.0
   ⋅   2.0    ⋅
  6.0   ⋅    8.0
   ⋅   3.0    ⋅
  8.0   ⋅   10.0
   ⋅   4.0    ⋅
 10.0   ⋅   12.0
   ⋅   5.0    ⋅


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
while indexing with slices will return a new `CommunityProfile`.
To get the values of a matrix slice, use `abundances` after indexing.

```jldoctest profiles
julia> comm[1,3]
4.0

julia> comm[6:8,3]
CommunityProfile{Float64, Taxon, MicrobiomeSample} with 3 features in 1 samples

Feature names:
g1, g2, g3

Sample names:
s3

julia> comm[1:3,3] |> abundances
3×1 SparseMatrixCSC{Float64, Int64} with 2 stored entries:
 4.0
  ⋅
 6.0
```

### Indexing with strings and regular expressions

It is often inconvenient to find the numerical index
of a particular feature or sample.
Instead, you can use strings or regular expressions
to get slices of a `CommunityProfile`,
which will match on the `String`  representation of the features or samples.
As with numerical indexing, if the index returns
a unique (feature, sample) pair,
it will return the abundance of that pair.
Otherwise, it will return a new CommunityProfile.

```jldoctest profiles
julia> comm["g__g1", "s1"]


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
3-element Vector{NamedTuple{(:sample, :subject)}}:
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
3-element Vector{NamedTuple{(:sample, :subject, :foo, :other)}}:
 (sample = "s1", subject = "kevin", foo = "bar", other = "Hello, World!")
 (sample = "s2", subject = "anika", foo = missing, other = "Goodbye!")
 (sample = "s3", subject = "annelle", foo = "baz", other = missing)
```

## Importing from tabular data

The `CommunityProfile` type is compatible with `Tables.jl`,
so should be able to import data from tabular sources
with the help of some other packages.

For example, using the [`DataFrames.jl`](https://github.com/JuliaData/DataFrames.jl) package
(which you can install at the REPL with `] add DataFrames`):

```jldoctest
julia> using DataFrames, Microbiome

julia> df = DataFrame(features = ["gene$i|s__Some_species" for i in 1:5], s1 = 1:2:10, s2 = 0:2:9)
5×3 DataFrame
 Row │ features               s1     s2
     │ String                 Int64  Int64
─────┼─────────────────────────────────────
   1 │ gene1|s__Some_species      1      0
   2 │ gene2|s__Some_species      3      2
   3 │ gene3|s__Some_species      5      4
   4 │ gene4|s__Some_species      7      6
   5 │ gene5|s__Some_species      9      8

julia> gfs = genefunction.(df.features)
5-element Vector{GeneFunction}:
 GeneFunction("gene1", Taxon("Some_species", :species))
 GeneFunction("gene2", Taxon("Some_species", :species))
 GeneFunction("gene3", Taxon("Some_species", :species))
 GeneFunction("gene4", Taxon("Some_species", :species))
 GeneFunction("gene5", Taxon("Some_species", :species))

julia> mss = MicrobiomeSample.(names(df)[2:end])
2-element Vector{MicrobiomeSample}:
 MicrobiomeSample("s1", {})
 MicrobiomeSample("s2", {})

julia> CommunityProfile(Matrix(df[!, 2:end]), gfs, mss)
CommunityProfile{Int64, GeneFunction, MicrobiomeSample} with 5 features in 2 samples

Feature names:
gene1, gene2, gene3, gene4, gene5

Sample names:
s1, s2
```

Alternatively, with a CSV file, and the [`CSV.jl`](https://github.com/JuliaData/CSV.jl) package:

```jldoctest csvexample
julia> println.(eachline("taxprof.csv"));
features,s1,s2
s__Escherichia_coli,1,0
s__Bifidobacterium_longum,3,2
s__Prevotella_copri,5,4

julia> using CSV, CSV.Tables

julia> tbl = CSV.read("taxprof.csv", Tables.columntable);

julia> txs = taxon.(tbl[1])
3-element Vector{Taxon}:
 Taxon("Escherichia_coli", :species)
 Taxon("Bifidobacterium_longum", :species)
 Taxon("Prevotella_copri", :species)

julia> mss = [MicrobiomeSample(string(k)) for k in keys(tbl)[2:end]]
2-element Vector{MicrobiomeSample}:
 MicrobiomeSample("s1", {})
 MicrobiomeSample("s2", {})

julia> mat =  hcat([tbl[i] for i in 2:length(tbl)]...)
3×2 Matrix{Int64}:
 1  0
 3  2
 5  4

julia> CommunityProfile(mat, txs, mss)
CommunityProfile{Int64, Taxon, MicrobiomeSample} with 3 features in 2 samples

Feature names:
Escherichia_coli, Bifidobacterium_longum, Prevotella_copri

Sample names:
s1, s2
```

You may also be interested in [`BiobakeryUtils.jl`](https://github.com/EcoJulia/BiobakeryUtils.jl)
which has convenient functions for reading in file types generated by
the `bioBakery` suite of computational tools.

## Types and Methods

```@docs
CommunityProfile
samples
features
samplenames
**featurenames**
commjoin
relativeabundance
relativeabundance!
present
prevalence
prevalence_filter
```
