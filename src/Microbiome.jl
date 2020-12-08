module Microbiome

# export
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

using EcoBase
using Statistics
using StatsBase

import Base: getindex, setindex, length

include("ecotranslations.jl")
include("abundances.jl")
include("distances.jl")

end  # module Microbiome
