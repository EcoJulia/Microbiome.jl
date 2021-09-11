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
## -- Convienience functions -- ##

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

profiletype(at::AbstractAbundanceTable) = eltype(features(at))
clades(at::AbstractAbundanceTable) = clade.(features(at))

Base.size(at::AbstractAbundanceTable, dims...) = size(at.aa, dims...)

# -- Indexing -- #

function Base.getindex(at::CommunityProfile, inds...)
    idx = at.aa[inds...]
    
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
    NamedTuple{(:features, Symbol.(samplenames(row))...)}((name(features(row)[1]), abundances(row)...))
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


## Metadata

"""
    metadata(cp::CommunityProfile)

Returns iterator of `NamedTuple` per sample, where keys are `:sample`
and each metadata key found in `cp`.
Samples without given metadata are filled with `missing`.

Returned values can be passed to any Tables.rowtable - compliant type,
eg `DataFrame`.
"""
function metadata(cp::CommunityProfile)
    ss = samples(cp)
    cols = unique(reduce(hcat, collect.(keys.(metadata.(samples(cp))))))
    return Tables.rowtable(merge((; sample=name(s)), 
                     NamedTuple(c => get(s, c, missing) for c in cols)
                    ) for s in ss)
end