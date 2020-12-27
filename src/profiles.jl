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

## -- Convienience functions -- ##

features(at::AbstractAbundanceTable) = axes(at.aa, 1) |> keys
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
