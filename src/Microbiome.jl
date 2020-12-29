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
       samples,
       clades,
       profiletype,
       featuretotals,
       sampletotals
   
# Abundances
export present,
       prevalence,
       relativeabundance!,
       relativeabundance       
#     filterabund

# Diversity
export ginisimpson,
       shannon,
       ginisimpson!,
       shannon!,
       present,
       prevalence,
       braycurtis,
       pcoa
       
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
include("abundances.jl")
include("distances.jl")

end  # module Microbiome
