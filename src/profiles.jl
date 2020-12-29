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
CommunityProfile{Float64, Taxon, MicrobiomeSample} with 10 things in 5 places

Thing names:
taxon1, taxon2, taxon3...taxon9, taxon10

Place names:
sample1, sample2, sample3, sample4, sample5

julia> comm["taxon1", "sample1"]
1.0

julia> comm[:,["sample1", "sample5"]]
CommunityProfile{Float64, Taxon, MicrobiomeSample} with 10 things in 2 places

Thing names:
taxon1, taxon2, taxon3...taxon9, taxon10

Place names:
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
        else
            error("invalid dimension name $dn")
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
# # todo
# EcoBase.thingoccurrences(at::AbstractAbundanceTable, things) = nothing
# EcoBase.placeoccurrences(at::AbstractAbundanceTable, places) = nothing

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

# ## -- Tables Interface -- ##

Tables.istable(::AbstractAbundanceTable) = true
Tables.columnaccess(::AbstractAbundanceTable) = true
Tables.rowaccess(::AbstractAbundanceTable) = true

Tables.getcolumn(at::AbstractAbundanceTable, i::Int) = i == 1 ? featurenames(at) : abundances(at[:, i-1])
Tables.getcolumn(at::AbstractAbundanceTable, i::Symbol) = i == :features ? featurenames(at) : abundances(at[:, string(i)])
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
