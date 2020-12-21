##################################
# Profiles - types and functions #
##################################

## -- Definitions and Constructors -- ##

abstract type MicrobiomeFeatures <: AbstractThings end
abstract type MicrobiomeSamples  <: AbstractPlaces{Nothing} end
# abstract type AbstractAbundanceTable{D, T, P} <: AbstractAssemblage{D, T, P} where {D <: Real, 
#                                                                             T <: MicrobiomeFeatures, 
#                                                                             P <: MicrobiomeSamples} end
abstract type AbstractAbundanceTable{D <: Real} end

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

mutable struct TaxonomicProfile{D} <: AbstractAbundanceTable{D}
    features::AbstractVector{Taxon}
    samples::AbstractVector{MicrobiomeSample}
    colindex::Dictionary{Symbol, Int}
    clades::AbstractVector{Union{Missing,Symbol}}
    abundances::SparseMatrixCSC{D}

    function TaxonomicProfile(tab, features, samples)
        D = eltype(tab)
        colindex = indexcols(samples)
        clades = clade.(features)
        return new{D}(features, samples, colindex, clades, tab)
    end
end


mutable struct FunctionalProfile{D} <: AbstractAbundanceTable{D}
    features::AbstractVector{GeneFunction}
    samples::AbstractVector{MicrobiomeSample}
    colindex::Dictionary{Symbol, Int}
    abundances::SparseMatrixCSC{D}
    
    function FunctionalProfile(tab, features, samples)
        D = eltype(tab)
        colindex = indexcols(samples)
        return new{D}(features, samples, colindex, tab)
    end
end

struct AbundanceTableRow{D} <: AbstractAbundanceTable{D}
    cols::NamedTuple

    function AbundanceTableRow(cols::NamedTuple)
        D = typeof(cols[2])
        return new{D}(cols)
    end
end

## -- Convienience functions -- ##

name(as::AbstractSample) = as.name
name(af::AbstractFeature) = af.name
clade(tax::Taxon) = tax.clade

taxon(gf::GeneFunction) = gf.taxon
hastaxon(gf::GeneFunction) = !ismissing(taxon(gf))

features(at::AbstractAbundanceTable) = at.features
samples(at::AbstractAbundanceTable) = at.samples
abundances(at::AbstractAbundanceTable) = at.abundances
abundances(atr::AbundanceTableRow) = values(atr.cols)[2:end]


Base.size(at::AbstractAbundanceTable, dims...) = size(abundances(at), dims...)
nthings(at::AbstractAbundanceTable) = size(at, 1)
nplaces(at::AbstractAbundanceTable) = size(at, 2)

# -- Indexing -- #

Base.getindex(at::AbstractAbundanceTable, ::Colon, i::Int) = i == 1 ? at.features : abundances(at)[:, i - 1]
Base.getindex(at::AbstractAbundanceTable, ::Colon, n::Symbol) = n == :features ? at.features : abundances(at)[:, at.colindex[n] - 1]
Base.getindex(atr::AbundanceTableRow, i::Union{Symbol, Int}) = atr.cols[i]

function Base.getindex(at::AbstractAbundanceTable, ri::Int, ::Colon)
    rowvals = abundances(at)[ri, :]
    return AbundanceTableRow((; :features => features(at)[ri], (Symbol(n) => rowvals[i] for (i,n) in enumerate(name(s) for s in at.samples))...))
end

function Base.getindex(at::AbstractAbundanceTable, f::AbstractString, ::Colon)
    rowidx = findfirst(r-> name(r) == f, features(at))
    return at[rowidx, :]
end

Base.getindex(at::AbstractAbundanceTable, row::AbstractString, col::Union{Symbol, Int}) = at[row][col]
Base.getindex(at::AbstractAbundanceTable, row::Int, col::Union{Symbol, Int}) = at[:, col][row]

## -- Views -- ##
## Basically lifted  all of this from SpatialEcology.jl,
## Thanks Michael Borregaard!

abstract type SubAbuncanceTable{D} <: AbstractAbundanceTable{D} end

mutable struct SubTaxonomicProfile{D} <: SubAbuncanceTable{D}
    features::AbstractVector{Taxon}
    samples::AbstractVector{MicrobiomeSample}
    colindex::Dictionary{Symbol, Int}
    clades::AbstractVector{Union{Missing,Symbol}}
    abundances::SparseMatrixCSC{D}
end

mutable struct SubFunctionalProfile{D} <: SubAbuncanceTable{D}
    features::AbstractVector{GeneFunction}
    samples::AbstractVector{MicrobiomeSample}
    colindex::Dictionary{Symbol, Int}
    abundances::SparseMatrixCSC{D}
end

## TODO
function view(abt::TaxonomicProfile; features = 1:nfeatures(com), samples = 1:nsamples(com))
    feat = EcoBase.asindices(sites, featurenames(abt))
    samp = EcoBase.asindices(samples, samplenames(abt))
    # SubComMatrix(view(abundances(atp), feat, samp), view(featurenames(abt), feat), view(samplenames(abt), samp)) #TODO change the order of these in the object to fit the array index order
end

## -- EcoBase Translations -- ##

EcoBase.thingnames(at::AbstractAbundanceTable) = name.(features(at))
EcoBase.placenames(at::AbstractAbundanceTable) = name.(samples(at))
EcoBase.occurrences(at::AbstractAbundanceTable) = abundances(at)

## -- Tables Interface -- ##

Tables.istable(::AbstractAbundanceTable) = true
Tables.columnaccess(::AbstractAbundanceTable) = true
Tables.rowaccess(::AbstractAbundanceTable) = true

Tables.getcolumn(at::AbstractAbundanceTable, i::Union{Int, Symbol}) = at[:, i]
Tables.getcolumn(atr::AbundanceTableRow, i::Union{Int, Symbol}) = atr[i]
Tables.columnnames(at::AbstractAbundanceTable) = [:features, Symbol.(name.(samples(at)))...]
Tables.columnnames(atr::AbundanceTableRow) = keys(atr.cols)

function Tables.schema(at::AbstractAbundanceTable)
    elt = eltype(abundances(at))
    coltypes = [eltype(features(at)), (elt for _ in 1:nsamples(at))...]
    return Tables.Schema(Tables.columnnames(at), coltypes)
end


Tables.columns(at::AbstractAbundanceTable) = (; (col => at[:, col] for col in Tables.columnnames(at))...)
Tables.rows(at::AbstractAbundanceTable) = [at[i, :] for i in 1:nfeatures(at)]
