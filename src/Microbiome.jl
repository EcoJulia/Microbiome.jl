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
       jaccard,
       hellinger,
       pcoa
       
using Statistics
using SparseArrays
using EcoBase
using AxisIndices
using Dictionaries
using NamedDims
using Tables
using Distances
using MultivariateStats

import Dictionaries: set!, unset!, insert!, delete!

include("ecobase.jl")
include("samples_features.jl")
include("profiles.jl")
include("diversity.jl")

end  # module Microbiome
