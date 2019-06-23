module Microbiome

export
    # Functions
    ## abundance
    abundancetable,
    filterabund,
    relativeabundance,
    relativeabundance!,
    rownormalize,
    rownormalize!,
    colnormalize,
    colnormalize!,
    nfeatures,
    getfeature,
    featurenames,
    featuretotals,
    nsamples,
    getsample,
    samplenames,
    sampletotals,
    ## Diversity
    ginisimpson,
    shannon,
    present,
    prevalence


using Reexport
@reexport using SpatialEcology

using Statistics
using StatsBase
using DataFrames

import SpatialEcology: @forward_func
import Base: getindex, setindex, length

include("ecotranslations.jl")
include("abundances.jl")
include("distances.jl")

end  # module Microbiome
