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

struct GeneFunction <: AbstractFeature
    name::String
    taxon::Union{Missing, Taxon}
end

function indexcols(fs::AbstractVector{<:AbstractSample})
    return Dictionary((:features, (Symbol(name(f)) for f in fs)...), 1:(length(fs) + 1))
end

mutable struct TaxonomicProfile <: AbundanceTable
    features::AbstractVector{Taxon}
    samples::AbstractVector{MicrobiomeSample}
    colindex::Dictionary{Symbol, Int}
    clades::AbstractVector{Union{Missing,Symbol}}
    abundances::SparseMatrixCSC

    function TaxonomicProfile(tab, features, samples)
        colindex = indexcols(samples)
        clades = clade.(features)
        return new(features, samples, colindex, clades, tab)
    end
end


mutable struct FunctionalProfile <: AbundanceTable
    features::AbstractVector{GeneFunction}
    samples::AbstractVector{MicrobiomeSample}
    colindex::Dictionary{Symbol, Int}
    abundances::SparseMatrixCSC
    function FunctionalProfile(tab, features, samples)
        colindex = indexcols(samples)
        new(features, samples, colindex, tab)
    end
end

struct AbundanceTableRow <: AbundanceTable
    cols::NamedTuple
end

## -- Convienience functions -- ##

name(as::AbstractSample) = as.name
name(af::AbstractFeature) = af.name
clade(tax::Taxon) = tax.clade

taxon(gf::GeneFunction) = gf.taxon
hastaxon(gf::GeneFunction) = !ismissing(taxon(gf))

features(at::AbundanceTable) = at.features
samples(at::AbundanceTable) = at.samples
abundances(at::AbundanceTable) = at.abundances
abundances(atr::AbundanceTableRow) = values(atr.cols)[2:end]


Base.size(at::AbundanceTable, dims...) = size(abundances(at), dims...)
nthings(at::AbundanceTable) = size(at, 1)
nplaces(at::AbundanceTable) = size(at, 2)

# -- Indexing -- #

Base.getindex(at::AbundanceTable, ::Colon, i::Int) = i == 1 ? at.features : abundances(at)[:, i - 1]
Base.getindex(at::AbundanceTable, ::Colon, n::Symbol) = n == :features ? at.features : abundances(at)[:, at.colindex[n] - 1]
Base.getindex(atr::AbundanceTableRow, i::Union{Symbol, Int}) = atr.cols[i]

function Base.getindex(at::AbundanceTable, ri::Int, ::Colon)
    rowvals = abundances(at)[ri, :]
    return AbundanceTableRow((; :features => features(at)[ri], (Symbol(n) => rowvals[i] for (i,n) in enumerate(name(s) for s in at.samples))...))
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
EcoBase.occurrences(at::AbundanceTable) = abundances(at)

## -- Tables Interface -- ##

Tables.istable(::AbundanceTable) = true
Tables.columnaccess(::AbundanceTable) = true
Tables.rowaccess(::AbundanceTable) = true

Tables.getcolumn(at::AbundanceTable, i::Union{Int, Symbol}) = at[:, i]
Tables.getcolumn(atr::AbundanceTableRow, i::Union{Int, Symbol}) = atr[i]
Tables.columnnames(at::AbundanceTable) = [:features, Symbol.(name.(samples(at)))...]
Tables.columnnames(atr::AbundanceTableRow) = keys(atr.cols)

function Tables.schema(at::AbundanceTable)
    elt = eltype(abundances(at))
    coltypes = [eltype(features(at)), (elt for _ in 1:nsamples(at))...]
    return Tables.Schema(Tables.columnnames(at), coltypes)
end


Tables.columns(at::AbundanceTable) = (; (col => at[:, col] for col in Tables.columnnames(at))...)
Tables.rows(at::AbundanceTable) = [at[i, :] for i in 1:nfeatures(at)]

# Tables.table(at::AbundanceTable) = 