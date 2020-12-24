module Microbiome

export nfeatures,
       featurenames,
       getfeature,
       nsamples,
       samplenames,
       getsample,
       AbstractFeature,
       AbstractSample,
       MicrobiomeSample,
       Taxon,
       GeneFunction,
       CommunityProfile,
       featuretype,
       name,
       clade,
       taxon,
       hastaxon,
       features,
       samples,
       abundances

#     # Functions
#     ## abundance
#     abundancetable,
#     filterabund,
#     relativeabundance,
#     relativeabundance!,
#     rownormalize,
#     rownormalize!,
#     colnormalize,
#     colnormalize!,
#     nfeatures,
#     getfeature,
#     featurenames,
#     featuretotals,
#     nsamples,
#     getsample,
#     samplenames,
#     sampletotals,
#     ## Diversity
#     ginisimpson,
#     shannon,
#     present,
#     prevalence

using Statistics
using StatsBase
using Tables
using Dictionaries
using SparseArrays
using EcoBase
import EcoBase: asindices
using Distances
using MultivariateStats
using AxisIndices
using NamedDims

import EcoBase: AbstractThings, AbstractPlaces, AbstractAssemblage,
                nthings, thingnames, thingoccurrences,
                nplaces, placenames, placeoccurrences

include("profiles.jl")
# include("tablesinterface.jl")
# include("abundances.jl")
include("distances.jl")

end  # module Microbiome
