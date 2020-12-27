const abundances = EcoBase.occurrences

const nfeatures = EcoBase.nthings
const featurenames = EcoBase.thingnames
const getfeature = EcoBase.thingoccurrences

const nsamples = EcoBase.nplaces
const samplenames = EcoBase.placenames
const getsample = EcoBase.placeoccurrences

abstract type AbstractFeature <: EcoBase.AbstractThings end
abstract type AbstractSample <: EcoBase.AbstractPlaces{Nothing} end