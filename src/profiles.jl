abstract type AbstractAbundanceTable{T <: Real, 
                                     F <: AbstractFeature, 
                                     S <: AbstractSample} <: EcoBase.AbstractAssemblage{T, F, S}
end

"""
    CommunityProfile{T, F, S} <: AbstractAbundanceTable{T, F, S}

An `AbstractAssemblage` from [EcoBase.jl](https://github.com/EcoJulia/EcoBase.jl)
that uses an `AxisArray` of a `SparseMatrixCSC` under the hood.

`CommunityProfile`s are tables with `AbstractFeature`-indexed rows and
`AbstractSample`-indexed columns.
Note - we can use the `name` of samples and features to index.

```jldoctest community
julia> txs = [Taxon("taxon\$i") for i in 1:10];

julia> mss = [MicrobiomeSample("sample\$i") for i in 1:5];

julia> mat = spzeros(10,5);

julia> for i in 1:5; mat[i,i] = 1.; end

julia> comm = CommunityProfile(mat, txs, mss)
CommunityProfile{Float64, Taxon, MicrobiomeSample} with 10 features in 5 samples

Feature names:
taxon1, taxon2, taxon3...taxon9, taxon10

Sample names:
sample1, sample2, sample3, sample4, sample5

julia> comm["taxon1", "sample1"]
1.0

julia> comm[:,["sample1", "sample5"]]
CommunityProfile{Float64, Taxon, MicrobiomeSample} with 10 features in 2 samples

Feature names:
taxon1, taxon2, taxon3...taxon9, taxon10

Sample names:
sample1, sample5

julia> comm[Taxon("taxon3", :kingdom), "sample1"]
0.0
```
"""
mutable struct CommunityProfile{T, F, S} <: AbstractAbundanceTable{T, F, S}
    aa::NamedAxisArray

    function CommunityProfile(aa::NamedAxisArray)
        @assert dimnames(aa) == (:features, :samples)
        T = eltype(parent(aa))
        F = eltype(keys(axes(aa, 1)))
        S = eltype(keys(axes(aa, 2)))
        return new{T, F, S}(aa)
    end
end

function CommunityProfile(tab::SparseMatrixCSC{<:Real}, 
                          features::AbstractVector{<:AbstractFeature},
                          samples::AbstractVector{<:AbstractSample})
    return CommunityProfile(NamedAxisArray(tab, features=features, samples=samples))
end

function CommunityProfile{T, F, S}(tab::SparseMatrixCSC{<:T},
                                   features::AbstractVector{F}, 
                                   samples::AbstractVector{S}) where {T, F, S}
    return CommunityProfile(tab, features, samples)
end

function CommunityProfile(tab::AbstractMatrix,
                          features::AbstractVector{<:AbstractFeature},
                          samples::AbstractVector{<:AbstractSample})
    return CommunityProfile(sparse(tab), features, samples)
end
## -- Convienience functions -- ##

function ==(p1::CommunityProfile, p2::CommunityProfile)
    return abundances(p1) == abundances(p2) && 
           samples(p1)    == samples(p2) &&
           features(p1)   == features(p2)
end

"""
    features(at::AbstractAbundanceTable)

Returns features in `at`. To get featurenames instead, use [`featurenames`](@ref).
"""
features(at::AbstractAbundanceTable) = axes(at.aa, 1) |> keys

"""
    samples(at::AbstractAbundanceTable)

Returns samples in `at`. To get samplenames instead, use [`samplenames`](@ref).
"""
samples(at::AbstractAbundanceTable) = axes(at.aa, 2) |> keys

"""
    samples(at::AbstractAbundanceTable, name::AbstractString)

Returns sample in `at` with name `name`.
"""
function samples(at::AbstractAbundanceTable, name::AbstractString)
    idx = findall(==(name), samplenames(at))
    length(idx) == 0 && throw(IndexError("No samples called $name"))
    length(idx) > 1 && throw(IndexError("More than one sample matches name $name"))
    return samples(at)[axes(at.aa, 2)][first(idx)]
end

profiletype(at::AbstractAbundanceTable) = eltype(features(at))
clades(at::AbstractAbundanceTable) = clade.(features(at))

Base.size(at::AbstractAbundanceTable, dims...) = size(at.aa, dims...)

Base.copy(at::AbstractAbundanceTable) = typeof(at)(copy(abundances(at)), copy(features(at)), deepcopy(samples(at)))

# -- Indexing -- #

function _index_profile(at, idx, inds)
    # single value - return that value
    ndims(idx) == 0 && return idx 
    # another table - return a new CommunityProfile with that table
    ndims(idx) == 2 && return CommunityProfile(idx)
    # a row or a column, figure out which, and make it 2D
    if ndims(idx) == 1
        dn = dimnames(idx)[1]
        # if it's a row...
        if dn == :samples
            return at[[inds[1]], inds[2]]
        # if it's a column
        elseif dn == :features
            return at[inds[1], [inds[2]]]
        end
    end
end

function _toinds(arr, inds::AbstractVector{<: Union{AbstractString, Regex}})
    return findall(a-> any(ind-> contains(a, ind), inds), arr)
end

# fall back ↑
_toinds(arr, ind::Union{AbstractString, Regex}) = _toinds(arr, [ind])

function _toinds(arr, inds::AbstractVector{<: Union{AbstractSample, AbstractFeature}})
    return findall(a-> any(==(a), inds), arr)
end

# fall back ↑
_toinds(arr, ind::Union{AbstractSample, AbstractFeature}) = _toinds(arr, [ind])

# if inds are integers, just return them
_toinds(_, ind::Int) = ind
_toinds(_, inds::AbstractVector{Int}) = inds

function Base.getindex(at::CommunityProfile, inds...)
    idx = at.aa[inds...]
    
    _index_profile(at, idx, inds)
end

function Base.getindex(at::CommunityProfile, rowind::Union{T, AbstractVector{<:T}} where T<:Union{AbstractString,Regex}, colind)
    rows = _toinds(featurenames(at), rowind)
    idx = at.aa[rows, colind]

    _index_profile(at, idx, (rows, colind))
end

function Base.getindex(at::CommunityProfile, rowind, colind::Union{T, AbstractVector{<:T}} where T<:Union{AbstractString,Regex})
    cols = _toinds(samplenames(at), colind)
    idx = at.aa[rowind, cols]

    _index_profile(at, idx, (rowind, cols))
end

function Base.getindex(at::CommunityProfile, rowind::Union{T, AbstractVector{<:T}} where T<:Union{AbstractString,Regex},
                                             colind::Union{S, AbstractVector{<:S}} where S<:Union{AbstractString,Regex})
    rows = _toinds(featurenames(at), rowind)
    at[rows, colind]
end

## -- EcoBase Translations -- ##
# see src/ecobase.jl for Microbiome function names
# thing => feature
# place => sample
# occurrences => abundances (or totals)

EcoBase.thingnames(at::AbstractAbundanceTable) = name.(features(at))
EcoBase.placenames(at::AbstractAbundanceTable) = name.(samples(at))
EcoBase.occurrences(at::AbstractAbundanceTable) = parent(parent(at.aa)) # first parent is the unnamed AxisArray
EcoBase.nthings(at::AbstractAbundanceTable) = size(at, 1)
EcoBase.nplaces(at::AbstractAbundanceTable) = size(at, 2)
## todo
# EcoBase.thingoccurrences(at::AbstractAbundanceTable, things) = nothing
# EcoBase.placeoccurrences(at::AbstractAbundanceTable, places) = nothing

# for custom printing
EcoBase.thingkind(asm::AbstractAbundanceTable) = "feature"
EcoBase.placekind(asm::AbstractAbundanceTable) = "sample"
## not needed for now
# EcoBase.thingkindplural(asm::AbstractAbundanceTable) = "$(thingkind(asm))s"
# EcoBase.placekindplural(asm::AbstractAbundanceTable) = "$(placekind(asm))s"

"""
    featuretotals(at::AbstractAbundanceTable)

Returns sum of each row (feature) in `at`.
Note, return value is a nfeatures x 1 `Matrix`, not a `Vector`.
If you need 1D `Vector`, use `vec(featuretotals(at))`.
"""
featuretotals(at::AbstractAbundanceTable) = sum(abundances(at), dims=2)

"""
    sampletotals(at::AbstractAbundanceTable)

Returns sum of each row (feature) in `at`.
Note, return value is a 1 x nsamples `Matrix`, not a `Vector`.
If you need 1D `Vector`, use `vec(sampletotals(at))`.
"""
sampletotals(at::AbstractAbundanceTable) = sum(abundances(at), dims=1)

## -- Tables Interface -- ##

Tables.istable(::AbstractAbundanceTable) = true
Tables.columnaccess(::AbstractAbundanceTable) = true
Tables.rowaccess(::AbstractAbundanceTable) = true

Tables.getcolumn(at::AbstractAbundanceTable, i::Int) = i == 1 ? featurenames(at) : abundances(at[:, i-1])
Tables.getcolumn(at::AbstractAbundanceTable, i::AbstractString) = i == "features" ? featurenames(at) : abundances(at[:, i])
Tables.getcolumn(at::AbstractAbundanceTable, i::Symbol) = Tables.getcolumn(at, string(i))

Tables.columnnames(at::AbstractAbundanceTable) = [:features, Symbol.(samplenames(at))...]

function Tables.schema(at::AbstractAbundanceTable)
    elt = eltype(abundances(at))
    coltypes = [eltype(features(at)), (elt for _ in 1:nsamples(at))...]
    return Tables.Schema(Tables.columnnames(at), coltypes)
end


Tables.columns(at::AbstractAbundanceTable) = (; (col => Tables.getcolumn(at, col) for col in Tables.columnnames(at))...)

function _makerow(row::AbstractAbundanceTable)
    size(row, 1) == 1 || error("Can't make row from table of size $(size(row))")
    NamedTuple{(:features, Symbol.(samplenames(row))...)}((first(featurenames(row)), abundances(row)...))
end

function _makerow(row::CommunityProfile{<:Real, <:GeneFunction, <:AbstractSample})
    size(row, 1) == 1 || error("Can't make row from table of size $(size(row))")
    gf = first(features(row))
    n = name(gf)
    hastaxon(gf) && (n *= string('|', name(taxon(gf))))
    NamedTuple{(:features, Symbol.(samplenames(row))...)}((n, abundances(row)...))
end


Tables.rows(at::AbstractAbundanceTable) = (_makerow(at[i, :]) for i in 1:nfeatures(at))

## -- Methods for absolute and relative abundances -- ##

# """
#     filterabund(abun::AbstractComMatrix, n::Int=minimum(10, nfeatures(abun)))

# Filter an abundance table to the top `n` features accross all samples

# This function also adds a row for "other", which sums the abundances of the
# remaining features.
# """
# function filterabund(abun::AbstractAbundanceTable, n::Int=minimum(10, nfeatures(abun)))
#     # TODO: add prevalence filter

#     totals = featuretotals(abun)

#     srt = sortperm(totals, rev=true)

#     newabun = getfeature(abun, srt[1:n])

#     remainder = [sum(occurrences(abun)[srt[n+1:end], i]) for i in 1:size(abun, 2)]'
#     newabun = vcat(newabun, remainder)
#     newrows = cat(featurenames(abun)[srt[1:n]], ["other"], dims=1)

#     return abundancetable(newabun, samplenames(abun), newrows)
# end

"""
    relativeabundance!(a::AbstractAbundanceTable; kind::Symbol=:fraction)

Normalize each sample in AbstractAbundanceTable to the sum of the sample.

By default, columns sum to 1.0.
Use `kind=:percent` for columns to sum to 100.
"""
function relativeabundance!(at::AbstractAbundanceTable; kind::Symbol=:fraction)
    in(kind, [:percent, :fraction]) || throw(ArgumentError("Invalid kind: $kind"))
    eltype(abundances(at)) <: AbstractFloat || throw(ArgumentError("relativeabundance! requires profile to have AbstractFloat eltype. Try relativeabundance instead"))
    abund = abundances(at)
    abund ./= sampletotals(at)
    kind == :percent && (abund .*= 100)
    dropzeros!(abund) # shouldn't need dropzeros - https://github.com/JuliaLang/julia/issues/39018
    return at
end

"""
    relativeabundance(at::AbstractAbundanceTable, kind::Symbol=:fraction)

Like [`relativeabundance!`](@ref), but does not mutate original.
"""
function relativeabundance(at::AbstractAbundanceTable, kind::Symbol=:fraction)
    comm = typeof(at)(float.(abundances(at)), deepcopy(features(at)), deepcopy(samples(at)))
    relativeabundance!(comm)
end

"""
    present(t::Union{Real, Missing}, minabundance::Real=0.0)
    present(at::AbstractAbundanceTable, minabundance::Real=0.0)

Check if a given (non-zero) value is greater than or equal to a minimum value.
If the minimum abundance is 0, just checks if value is non-zero.

If used on an `AbstractAbundanceTable`, returns a sparse boolean matrix of the same size.
"""
function present(t::Real, minabundance::Real=0.0)
    (minabundance >= 0 && t >= 0) || throw(DomainError("Only defined for positive values"))
    t == 0 ? false : t >= minabundance
end

present(::Missing, m) = missing

function present(at::AbstractAbundanceTable, minabundance::Real=0.0)
    mat = spzeros(Bool, Ssize(at)...)
    for i in eachindex(mat)
        mat[i] = present(at[i], minabundance)
    end
    return mat
end


"""
    prevalence(a::AbstractArray{<:Real}, minabundance::Real=0.0)
    prevalence(at::AbstractAbundanceTable, minabundance::Real=0.0)

Return the fraction of values that are greater than or equal to a minimum.
If the minimum abundance is 0, returns the fraction of non-zero values.

If used on an `AbstractAbundanceTable`,
returns a prevalence value for each `feature` accross the `sample`s.
"""
prevalence(a::AbstractArray{<:Real}, minabundance::Real=0.0) = mean(x-> present(x, minabundance), a)

# makes it work for any iterable
prevalence(a, minabundance::Real=0.0) = mean(x-> present(x, minabundance), (y for y in a))

function prevalence(at::AbstractAbundanceTable, minabundance::Real=0.0)
    mean(x-> present(x, minabundance), abundances(at), dims=2)
end

"""
    prevalence_filter(comm::AbstractAbundanceTable; minabundance=0.0; minprevalence=0.05, renorm=false)

Return a filtered `CommunityProfile` where features with prevalence lower than `minprevalence` are removed.
By default, a feature is considered "present" if > 0, but this can be changed by setting `minabundance`.

Optionally, set `renorm = true` to calculate relative abundances after low prevalence features are removed.

```jldoctest
julia> comm = CommunityProfile(sparse([3 0 1 # 0.33, assuming minabundance 2
                                       2 2 2 # 1.0
                                       0 0 1 # 0.0
                                       2 0 0 # 0.33
                                       ]),
                               [Taxon(string(i)) for i in 1:4],
                               [MicrobiomeSample(string(i)) for i in 1:3]);

julia> prevalence_filter(comm, minabundance=2, minprevalence=0.3) 
CommunityProfile{Int64, Taxon, MicrobiomeSample} with 3 features in 3 samples

Feature names:
1, 2, 4

Sample names:
1, 2, 3

julia> prevalence_filter(comm, minabundance=2, minprevalence=0.4)
CommunityProfile{Int64, Taxon, MicrobiomeSample} with 1 features in 3 samples

Feature names:
2

Sample names:
1, 2, 3

julia> prevalence_filter(comm, minabundance=3, minprevalence=0.3)
CommunityProfile{Int64, Taxon, MicrobiomeSample} with 1 features in 3 samples

Feature names:
1

Sample names:
1, 2, 3
```
"""
function prevalence_filter(comm::AbstractAbundanceTable; minabundance=0.0, minprevalence=0.05, renorm=false)
    comm = comm[vec(prevalence(comm, minabundance) .>= minprevalence), :]
    return renorm ? relativeabundance(comm) : comm
end

"""
    cladefilter(comm::AbstractAbundanceTable, cl::Union{Symbol, Int}; keepempty=false)

Return a copy of `comm`, where only rows that have `clade(feature) == cl` are kept.
Use `keepempty = true` to also keep features that don't have a `clade` (eg "UNIDENTIFIED").

Examples
≡≡≡≡≡≡≡≡≡≡

```jldoctest
julia> features(comm)
10-element Vector{Taxon}:
 Taxon("taxon1", :domain)
 Taxon("taxon2", :kingdom)
 Taxon("taxon3", :phylum)
 Taxon("taxon4", :class)
 Taxon("taxon5", :order)
 Taxon("taxon6", :family)
 Taxon("taxon7", :genus)
 Taxon("taxon8", :species)
 Taxon("taxon9", :subspecies)
 Taxon("taxon10", missing)

julia> features(cladefilter(comm, :species))
 1-element Vector{Taxon}:
  Taxon("taxon8", :species)

julia> features(cladefilter(comm, :genus; keepempty = true))
  2-element Vector{Taxon}:
   Taxon("taxon7", :genus)
   Taxon("taxon10", missing)
```
"""
function cladefilter(comm::AbstractAbundanceTable, cl::Symbol; keepempty=false)
    in(cl, keys(_clades)) ||  error("Invalid clade $cl, must be one of $(keys(_clades))")
    if keepempty
        return filter(f-> !hasclade(f) || clade(f) == cl, comm)
    else
        return filter(f-> hasclade(f) && clade(f) == cl, comm)
    end
end

function cladefilter(comm::AbstractAbundanceTable, clade::Int; keepempty=false)
    0 <= clade <= 9 ||  error("Invalid clade $clade, must be one of $_clades")
    return cladefilter(comm, keys(_clades)[clade+1]; keepempty)
end


## Metadata

"""
    metadata(commp::CommunityProfile)

Returns iterator of `NamedTuple` per sample, where keys are `:sample`
and each metadata key found in `commp`.
Samples without given metadata are filled with `missing`.

Returned values can be passed to any Tables.rowtable - compliant type,
eg `DataFrame`.
```
"""
function metadata(commp::CommunityProfile)
    ss = samples(commp)
    cols = unique(reduce(vcat, collect.(keys.(metadata.(samples(commp))))))
    return Tables.rowtable(merge((; sample=name(s)), 
                     NamedTuple(c => get(s, c, missing) for c in cols)
                    ) for s in ss)
end

"""
    add_metadata!(commp::CommunityProfile, samplename::AbstractString, md::Union{AbstractDict,NamedTuple}; overwrite=false)

Add metadata (in the form of an `AbstractDict` or `NamedTuple`) to the `MicrobiomeSample` in `commp` with name `samplename`.
For `AbstractDict`s, all keys must be `Symbol`s. 

The function will fail if any of the keys in `md` already exist in the `MicrobiomeSample`,
unless `overwrite=true` is used.

Examples
≡≡≡≡≡≡≡≡≡≡

```jldoctest
julia> metadata(comm)
3-element Vector{NamedTuple{(:sample,), Tuple{String}}}:
 (sample = "sample1",)
 (sample = "sample2",)
 (sample = "sample3",)

julia> add_metadata!(comm, "sample1", Dict(:subjectname=>"kevin", :age=>37))

julia> metadata(comm)
3-element Vector{NamedTuple{(:sample, :subjectname, :age), T} where T<:Tuple}:
 (sample = "sample1", subjectname = "kevin", age = 37)
 (sample = "sample2", subjectname = missing, age = missing)
 (sample = "sample3", subjectname = missing, age = missing)
"""
function add_metadata!(commp::CommunityProfile, samplename::AbstractString, md::Union{AbstractDict,NamedTuple}; overwrite=false)
    s = samples(commp, samplename)
    if !overwrite
        length(keys(md) ∩ keys(metadata(s))) == 0 || throw(IndexError("Adding this metadata would overwrite existing values. Use `overwrite=true` to proceed anyway"))
    end
    
    for key in keys(md)
        value = md[key]
        overwrite ? set!(s, key, value) : insert!(s, key, value)
    end
    return nothing
end

"""
    add_metadata!(commp::CommunityProfile, samplecol::Symbol, md; overwrite=false)

Add metadata (in the form of a `Tables.jl` table) a `CommunityProfile`.
One column (`samplecol`) should contain sample names that exist in `commp`,
and other columns should contain metadata that will be added to the metadata of each sample.

The function will fail if any of the column names in `md` already exist as metadata in any of the `MicrobiomeSample`s,
unless `overwrite=true` is used.

Examples
≡≡≡≡≡≡≡≡≡≡

```jldoctest
julia> metadata(comm)
3-element Vector{NamedTuple{(:sample,), Tuple{String}}}:
 (sample = "sample1",)
 (sample = "sample2",)
 (sample = "sample3",)

julia> md_table = [(id="sample1", something=5,  newthing="bar"),
                   (id="sample2", something=10, newthing="baz"),
                   (id="sample3", something=42, newthing="fuz")];

julia> add_metadata!(comm, :id, md_table)

julia> metadata(comm)
3-element Vector{NamedTuple{(:sample, :something, :newthing), Tuple{String, Int6
4, String}}}:
 (sample = "sample1", something = 5, newthing = "bar")
 (sample = "sample2", something = 10, newthing = "baz")
 (sample = "sample3", something = 42, newthing = "fuz")
 ```
"""
function add_metadata!(commp::CommunityProfile, samplecol::Symbol, md; overwrite = false)
    Tables.istable(md) || throw(ArgumentError("Metadata must be a Tables.table"))
    for row in Tables.rows(md)
        row[samplecol] in samplenames(commp) || throw(IndexError("Sample '$(row[samplecol])' not found in CommunityProfile"))
        sample = samples(commp, row[samplecol])
        if !overwrite
            any(k-> haskey(sample, k), Tables.columnnames(row)) && throw(IndexError("Adding this metadata would overwrite existing values. Use `overwrite=true` to proceed anyway"))
        end
    end

    for row in Tables.rows(md)
        rowmd = Dict(col=> row[col] for col in Tables.columnnames(row) if col != samplecol)
        add_metadata!(commp, row[samplecol], rowmd; overwrite)
    end
    return nothing
end

"""
    set!(commp::CommunityProfile, sample::AbstractString, prop::Symbol, val)

Update or insert a value `val` to the metadata of `sample` in the CommunityProfile `commp` using a Symbol `prop`. 
If you want an error to be thrown if the value already exists, use [`insert!`](@ref).

Examples
≡≡≡≡≡≡≡≡≡≡

```jldoctest
julia> set!(comm, "sample1", :something, 1.0)

julia> first(metadata(comm))[:something]
1.0
 ```
"""
function set!(commp::CommunityProfile, sample::AbstractString, prop::Symbol, val)
    sample = samples(commp, sample)
    prop in _restricted_fields(sample) && error("Cannot set! $prop for $(typeof(sample)).")
    set!(sample.metadata, prop, val)
    return sample
end

"""
    unset!(commp::CommunityProfile, sample::AbstractString, prop::Symbol)

Delete a metadata entry in `sample` from CommunityProfile `commp` using the Symbol `prop`. 
If you want an error to be thrown if the value does not exist, use [`delete!`](@ref).

Examples
≡≡≡≡≡≡≡≡≡≡

```jldoctest
julia> unset!(comm, "sample1", :something) 

julia> !haskey(comm, "sample1", :something)
true
 ```
"""
function unset!(commp::CommunityProfile, sample::AbstractString, prop::Symbol)
    sample = samples(commp, sample)
    prop in _restricted_fields(sample) && error("Cannot unset! $prop for $(typeof(sample)).")
    unset!(sample.metadata, prop)
    return sample
end

"""
    insert!(commp::CommunityProfile, sample::AbstractString, prop::Symbol, val)

Insert a value `val` to the metadata of `sample` in a CommunityProfile `commp` using a Symbol `prop`, 
and it will throw an error if `prop` exists. 
If you don't want an error to be thrown if the value exists, use [`set!`](@ref).


Examples
≡≡≡≡≡≡≡≡≡≡

```jldoctest
julia> insert!(comm, "sample1", :something, 3.0) 

julia> get(comm, "sample1", :something, 3.0)
3.0
 ```
"""
function insert!(commp::CommunityProfile, sample::AbstractString, prop::Symbol, val)
    sample = samples(commp, sample)
    prop in _restricted_fields(sample) && error("Cannot insert! $prop for $(typeof(sample)).")
    insert!(sample.metadata, prop, val)
    return sample
end

"""
    delete!(commp::CommunityProfile, sample::AbstractString, prop::Symbol)

Delete a metadata entry in `sample` from CommunityProfile `commp` using the Symbol `prop` if it exists, or throw an error otherwise.
If you don't want an error to be thrown if the value does not exist, use [`unset!`](@ref).

Examples
≡≡≡≡≡≡≡≡≡≡

```jldoctest
julia> delete!(comm, "sample1", :something) 

julia> !haskey(comm, "sample1", :something)
true
 ```
"""
function delete!(commp::CommunityProfile, sample::AbstractString, prop::Symbol)
    sample = samples(commp, sample)
    prop in _restricted_fields(sample) && error("Cannot delete! $prop for $(typeof(sample)).")
    delete!(sample.metadata, prop)
    return sample
end

"""
    keys(commp::CommunityProfile, sample::AbstractString)

Return an iterator over all keys of the metadata attached to `sample` in a CommunityProfile `commp`. 
`collect(keys(commp, sample))` returns an array of keys. 

Examples
≡≡≡≡≡≡≡≡≡≡

```jldoctest
julia> add_metadata!(comm, "sample1", Dict(:subjectname=>"kevin", :age=>37))

julia> collect(keys(comm, "sample1"))
2-element Vector{Symbol}:
 :subjectname
 :age
```
"""
Base.keys(commp::CommunityProfile, sample::AbstractString) = keys(metadata(samples(commp, sample)))

"""
    haskey(commp::CommunityProfile, sample::AbstractString, key::Symbol)

Determine whether the metadata of `sample` in a CommunityProfile `commp` has a mapping for a given `key`. 
Use `!haskey` to determine whether a `sample` in a CommunityProfile doesn't have a mapping for a given `key`

Examples
≡≡≡≡≡≡≡≡≡≡

```jldoctest
julia> set!(comm, "sample1", :something, 1.0)

julia> haskey(comm, "sample1", :something)
true

julia> delete!(comm, "sample1", :something) 

julia> !haskey(comm, "sample1", :something)
true
 ```
"""
Base.haskey(commp::CommunityProfile, sample::AbstractString, key::Symbol) = in(key, keys(samples(commp, sample)))

"""
    get(commp::CommunityProfile, sample::AbstractString, key::Symbol, default)

Return the value of the metadata in a `sample` stored for the given `key`, or the given `default` value if no mapping for the key is present.

Examples
≡≡≡≡≡≡≡≡≡≡

```jldoctest
julia> get(comm, "sample1", :something, 42)
42 

julia> insert!(comm, "sample1", :something, 3.0) 

julia> get(comm, "sample1", :something, 42)
3.0
 ```
"""
Base.get(commp::CommunityProfile, sample::AbstractString, key::Symbol, default) = get(metadata(samples(commp, sample)), key, default)


"""
    filter(f, comm::CommunityProfile)

Apply `f` to the features of `comm`,
and return a copy where `f(feature)` is `true`.

Examples
≡≡≡≡≡≡≡≡≡≡

```jldoctest
julia> features(comm)
3-element Vector{GeneFunction}:
 GeneFunction("gene1", Taxon("tax1", :species))
 GeneFunction("gene1", Taxon("tax2", :genus))
 GeneFunction("gene2", missing)

julia> features(filter(hastaxon, comm))
2-element Vector{GeneFunction}:
 GeneFunction("gene1", Taxon("tax1", :species))
 GeneFunction("gene1", Taxon("tax2", :genus))
```
"""
function Base.filter(f::Function, commp::CommunityProfile)
    ridx = findall(f, features(commp))
    isempty(ridx) && error("Can't return empty profile")
    return copy(commp[ridx, :])
end