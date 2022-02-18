module Microbiome

# Samples and Features
export MicrobiomeSample,
       Taxon,
       GeneFunction,
       Metabolite,
       metadata,
       set!,
       unset!,
       insert!,
       delete!,
       name,
       taxrank,
       hasrank,
       taxon,
       hastaxon,
       genefunction,
       commonname,
       masscharge,
       retentiontime

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
       ranks,
       rankfilter,
       profiletype,
       featuretotals,
       sampletotals,
       commjoin,
       metadata
   
# Abundances
export present,
       prevalence,
       relativeabundance!,
       relativeabundance,
       prevalence_filter
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
using Dictionaries
using Tables
using Distances
using MultivariateStats
using ReTest

import Dictionaries: set!, unset!, insert!, delete!
import Base: ==

include("ecobase.jl")
include("samples.jl")
include("features.jl")
include("profiles.jl")
include("diversity.jl")
include("comm_joins.jl")

end  # module Microbiome
