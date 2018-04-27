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
    ## similarity
    getdm,
    getrowdm,
    pcoa,
    eigenvalue,
    principalcoord,
    variance,
    optimalorder,
    optimalorder!,
    ## plotting
    hclustplot,
    annotationbar,
    ## utils
    metaphlan_import,
    panphlan_calcs,
    bysample,
    taxfilter,
    taxfilter!

using Reexport
@reexport using SpatialEcology
@reexport using Distances

using RecipesBase
using StatPlots
using StatsBase
using Colors
using DataFrames
using FileIO
using CSVFiles
using MicroLogging

import SpatialEcology.@forward_func
import SpatialEcology.summary
import SpatialEcology.show
import Clustering: Hclust, hclust
import Base: getindex, setindex, length

include("ecotranslations.jl")
include("abundances.jl")
include("distances.jl")
include("leafordering.jl")
include("plotting.jl")
include("biobakery_utils.jl")


end  # module Microbiome
