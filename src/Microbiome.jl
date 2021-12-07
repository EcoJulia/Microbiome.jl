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
       taxrank,
       hasrank,
       taxon,
       hastaxon,
       genefunction

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
using AxisIndices
using Dictionaries
using NamedDims
using Tables
using Distances
using MultivariateStats
using ReTest

import Dictionaries: set!, unset!, insert!, delete!
import Base: ==

@testset "test" begin
    @test true
end

include("ecobase.jl")
include("samples.jl")
include("features.jl")
include("profiles.jl")
include("diversity.jl")
include("comm_joins.jl")

end  # module Microbiome

"Microbiome.CommunityProfile{Float64, Microbiome.Taxon, Microbiome.MicrobiomeSample} with 10 features in 5 samples\n\nFeature names:\ntaxon1, taxon2, taxon3...taxon9, taxon10\n\nSample names:\nsample1, sample2, sample3, sample4, sample5\n\n" == 
"CommunityProfile{Float64, Taxon, MicrobiomeSample} with 10 features in 5 samples\n\nFeature names:\ntaxon1, taxon2, taxon3...taxon9, taxon10\n\nSample names:\nsample1, sample2, sample3, sample4, sample5\n\n"