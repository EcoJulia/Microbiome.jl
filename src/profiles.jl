##################################
# Profiles - types and functions #
##################################

## -- Definitions and Constructors -- ##

abstract type MicrobiomeFeatures <: AbstractThings end
abstract type MicrobiomeSamples  <: AbstractPlaces{Nothing} end
# abstract type AbundanceTable{D, T, P} <: AbstractAssemblage{D, T, P} where {D <: Real, 
#                                                                             T <: MicrobiomeFeatures, 
#                                                                             P <: MicrobiomeSamples} end
abstract type AbundanceTable end

const nfeatures = nthings
const featurenames = thingnames
const getfeature = thingoccurrences
# const featuretotals = speciestotals

const nsamples = nplaces
const samplenames = placenames
const getsample = placeoccurrences
# const sampletotals = sitetotals

abstract type AbstractFeature end
abstract type AbstractSample end

struct MicrobiomeSample
    name::String
    metadata::Dictionary{Symbol, T} where {T}
end

# Base.convert(MicrobiomeSample, s::String) = MicrobiomeSample(s, Dictionary{Symbol, Any}())

struct Taxon <: AbstractFeature
    name::String
    taxonlevel::Union{Missing, Symbol}
end


struct GeneFunction <: AbstractFeature
    name::String
    taxon::Union{Missing, Taxon}
end

function indexcols(fs::AbstractVector{<:AbstractFeature})
    return Dictionary((:features, (Symbol(name(f)) for f in fs)...), 1:(length(fs) + 1))
end

mutable struct TaxonomicProfile <: AbundanceTable
    features::AbstractVector{Taxon}
    samples::AbstractVector{MicrobiomeSample}
    colindex::Dictionary{Symbol, Int}
    taxonlevels::AbstractVector{Union{Missing,Symbol}}
    abundances::SparseMatrixCSC

    function TaxonomicProfile(tab, features, samples)
        colindex = indexcols(features)
        taxonlevels = taxonlevels.(features)
    end
end


mutable struct FunctionalProfile
    features::AbstractVector{GeneFunction}
    samples::AbstractVector{MicrobiomeSample}
    colindex::Dictionary{Symbol, Int}
    abundances::SparseMatrixCSC
end

struct AbundanceTableRow <: AbundanceTable
    cols::NamedTuple
end

## -- Convienience functions -- ##

name(af::AbstractFeature) = af.name
level(tax::Taxon) = tax.taxonlevel

taxon(gf::GeneFunction) = gf.taxon
hastaxon(gf::GeneFunction) = !ismissing(genetaxon(gf))

features(at::AbundanceTable) = at.features
samples(at::AbundanceTable) = at.samples
abundancetable(at::AbundanceTable) = at.abundances

Base.size(at::AbundanceTable, dims...) = size(abundancetable(at), dims...)
nthings(at::AbundanceTable) = size(at, 1)
nplaces(at::AbundanceTable) = size(at, 2)

# -- Indexing -- #

Base.getindex(at::AbundanceTable, ::Colon, i::Int) = n == 1 ? at.features : abundancetable(at)[:, i - 1]
Base.getindex(at::AbundanceTable, ::Colon, n::Symbol) = n == :features ? at.features : abundancetable(at)[:, at.colmap[n]]
Base.getindex(atr::AbundanceTableRow, i::Union{Symbol, Int}) = atr.cols[i]

function Base.getindex(at::AbundanceTable, ri::Int, ::Colon)
    rowvals = abundancetable(at)[ri, :]
    return AbundanceTableRow(; :features => f, (Symbol(n) => rowvals[i] for (i,n) in enumerate(name(s) for s in at.samples))...)
end

function Base.getindex(at::AbundanceTable, f::AbstractString, ::Colon)
    rowidx = findfirst(r-> name(r) == f, features(at))
    return at[rowidx, :]
end

Base.getindex(at::AbundanceTable, row::AbstractString, col::Union{Symbol, Int}) = at[row][col]
Base.getindex(at::AbundanceTable, row::Int, col::Union{Symbol, Int}) = at[:, col][row]

## -- EcoBase Translations -- ##

EcoBase.thingnames(at::AbundanceTable) = name.(features(at))
EcoBase.placenames(at::AbundanceTable) = name.(samples(at))
EcoBase.occurrences(at::AbundanceTable) = abundancetable(at)

## -- Tables Interface -- ##

Tables.istable(::AbundanceTable) = true
Tables.columnaccess(::AbundanceTable) = true
Tables.rowaccess(::AbundanceTable) = true

Tables.getcolumn(at::AbundanceTable, i::Union{Int, Symbol}) = at[:, i]
Tables.getcolumn(atr::AbundanceTableRow, i::Union{Int, Symbol}) = atr[i]
Tables.columnnames(at::AbundanceTable) = [:features, Symbol.(name.(samples(at)))...]
Tables.columnnames(atr::AbundanceTableRow) = keys(atr.cols)

Tables.columns(at::AbundanceTable) = (at[:, i] for i in 1:(nfeatures(at) + 1))
Tables.rows(at::AbundanceTable) = (at[i, :] for i in 1:nfeatures(at))