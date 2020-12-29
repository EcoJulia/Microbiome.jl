"""
    abundances(at::AbstractAbundanceTable)

Get the underlying sparse matrix of an `AbstractAbundanceTable`.
Note that this does not copy - any modifications to this matrix
will update the parent.
"""
const abundances = EcoBase.occurrences

const nfeatures = EcoBase.nthings

"""
    featurenames(at::AbstractAbundanceTable)

Get a vector of feature names from `at`, equivalent to `name.(features(at))`
"""
const featurenames = EcoBase.thingnames
const getfeature = EcoBase.thingoccurrences

const nsamples = EcoBase.nplaces

"""
    samplenames(at::AbstractAbundanceTable)

Get a vector of sample names from `at`, equivalent to `name.(samples(at))`
"""
const samplenames = EcoBase.placenames
const getsample = EcoBase.placeoccurrences

abstract type AbstractFeature <: EcoBase.AbstractThings end
abstract type AbstractSample <: EcoBase.AbstractPlaces{Nothing} end