module Microbiome

# Samples and Features
export MicrobiomeSample,
       Taxon,
       GeneFunction,
       metadata,
       set!,
       unset!,
       insert!,
       delete!,
       name,
       clade,
       hasclade,
       taxon,
       hastaxon

# EcoBase Translations
export abundances,
       nfeatures,
       featurenames,
       getfeature,
       nsamples,
       samplenames,
       getsample,
       AbstractFeature,
       AbstractSample
    #    featureabundances,
    #    sampleabundances

# Profiles
export CommunityProfile,
       featuretype,
       features,
       samples
    #    featuretotals,
    #    sampletotals
   
       

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
using Distances
using MultivariateStats
using AxisIndices
using NamedDims
using Dictionaries

import Dictionaries: set!, unset!, insert!, delete!

include("ecobase.jl")
include("samples_features.jl")
include("profiles.jl")
# include("tablesinterface.jl")
# include("abundances.jl")
include("distances.jl")

end  # module Microbiome
