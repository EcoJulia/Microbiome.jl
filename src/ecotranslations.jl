# # Interface with EcoBase

abstract type MicrobiomeFeatures <: AbstractThings end
abstract type MicrobiomeSamples  <: AbstractPlaces end
abstract type AbundanceTable <: AbstractAssemblage end

const nfeatures = nthings
const featurenames = thingnames
const getfeature = thingoccurrences
# const featuretotals = speciestotals

const nsamples = nplaces
const samplenames = placenames
const getsample = placeoccurrences
# const sampletotals = sitetotals

abstract type AbstractFeature end
name(af::AbstractFeature) = af.name

struct Taxon <: AbstractFeature
    name::String
    taxonlevel::Symbol
end

level(tax::Taxon) = tax.taxonlevel

struct GeneFunction <: AbstractFeature
    name::String
    taxon::Union{Missing, Taxon}
end

genetaxon(gf::GeneFunction) = gf.taxon
taxonname(gf::GeneFunction) = name(genetaxon(gf))
hastaxon(gf::GeneFunction) = !ismissing(genetaxon(gf))

mutable struct TaxonomicProfile{T}
    features::AbstractVector{Taxon}
    samplenames::AbstractVector{String}
    taxonlevels::AbstractVector{Union{Missing,Symbol}}
    table::SparseMatrixCSC{T}
end

mutable struct FunctionalProfile{T}
    features::AbstractVector{GeneFunction}
    samplenames::AbstractVector{String}
    table::SparseMatrixCSC{T}
end
