##################################
# Profiles - types and functions #
##################################

## -- Definitions and Constructors -- ##

const abundances = EcoBase.occurrences

const nfeatures = EcoBase.nthings
const featurenames = EcoBase.thingnames
const getfeature = EcoBase.thingoccurrences
# const featuretotals = speciestotals

const nsamples = EcoBase.nplaces
const samplenames = EcoBase.placenames
const getsample = EcoBase.placeoccurrences
# const sampletotals = sitetotals

abstract type AbstractFeature <: AbstractThings end
abstract type AbstractSample <: AbstractPlaces{Nothing} end

featuretype(af::AbstractFeature) = error("No feature type defined for $(typeof(af))")

abstract type AbstractAbundanceTable{T <: Real, 
                                     F <: AbstractFeature, 
                                     S <: AbstractSample} <: EcoBase.AbstractAssemblage{T, F, S}
end


struct MicrobiomeSample <: AbstractSample
    name::String
    metadata::Dictionary{Symbol, T} where {T <: Any} # currently non-functional
end

MicrobiomeSample(n::AbstractString) = MicrobiomeSample(n, Dictionary{Symbol, Any}())
# Base.convert(::Type{MicrobiomeSample}, s::AbstractString) = MicrobiomeSample(s)

const _clades = (
    domain     = 0,
    kingdom    = 1,
    phylum     = 2,
    class      = 3,
    order      = 4,
    family     = 5,
    genus      = 6,
    species    = 7,
    subspecies = 8,
    strain     = 9
)

struct Taxon <: AbstractFeature
    name::String
    clade::Union{Missing, Symbol}
    
    Taxon(s::AbstractString, ::Missing) = new(s, missing)
    Taxon(s::AbstractString, clade::Symbol) = in(clade, keys(_clades)) ? 
                                                        new(s, clade)  :
                                                        error("Invalid clade $clade, must be one of $(keys(_clades))")
end

Taxon(n::AbstractString, clade::Int) = 0 <= clade <= 9 ?
                                            Taxon(n, keys(_clades)[clade+1]) :
                                            error("Invalid clade $clade, must be one of $_clades")
Taxon(n::AbstractString) = Taxon(n, missing)
featuretype(::Taxon) = :taxon

struct GeneFunction <: AbstractFeature
    name::String
    taxon::Union{Missing, Taxon}
end

GeneFunction(n::AbstractString) = GeneFunction(n, missing)
featuretype(::GeneFunction) = :genefunction

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

featuretype(::AbstractAbundanceTable{T, F, S}) where {T, F <: Taxon, S} = :taxon
featuretype(::AbstractAbundanceTable{T, F, S}) where {T, F <: GeneFunction, S} = :genefunction

## -- Convienience functions -- ##

name(as::AbstractSample) = as.name
name(af::AbstractFeature) = af.name
Base.:(==)(s1::T, s2::T) where {T <: Union{AbstractSample, AbstractSample}} = name(s1) == name(s2)

clade(::Missing) = missing
clade(tax::Taxon) = tax.clade

taxon(gf::GeneFunction) = gf.taxon
clade(gf::GeneFunction) = clade(taxon(gf))
hastaxon(gf::GeneFunction) = !ismissing(taxon(gf))
hasclade(af::AbstractFeature) = !ismissing(clade(af))

features(at::AbstractAbundanceTable) = axes(at.aa, 1) |> keys
samples(at::AbstractAbundanceTable) = axes(at.aa, 2) |> keys

profiletype(at::AbstractAbundanceTable) = eltype(features(at))
clades(at::AbstractAbundanceTable) = clade.(features(at))


Base.size(at::AbstractAbundanceTable, dims...) = size(at.aa, dims...)
nthings(at::AbstractAbundanceTable) = size(at, 1)
nplaces(at::AbstractAbundanceTable) = size(at, 2)

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

EcoBase.thingnames(at::AbstractAbundanceTable) = name.(features(at))
EcoBase.placenames(at::AbstractAbundanceTable) = name.(samples(at))
EcoBase.occurrences(at::AbstractAbundanceTable) = parent(parent(at.aa)) # first parent is the unnamed AxisArray

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
