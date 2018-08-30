module Microbiome

export
    # Types
    DistanceMatrix,
    PCoA,

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
    ## distances
    getdm,
    getrowdm,
    pcoa,
    eigenvalue,
    principalcoord,
    variance,
    optimalorder,
    optimalorder!,
    ginisimpson,
    shannon,
    present,
    prevalence


using Reexport
@reexport using SpatialEcology
@reexport using Distances

using LinearAlgebra
using Statistics
using StatsBase
using DataFrames
using Clustering

import SpatialEcology: @forward_func
import Base: getindex, setindex, length

include("ecotranslations.jl")
include("abundances.jl")
include("distances.jl")
include("leafordering.jl")

end  # module Microbiome
