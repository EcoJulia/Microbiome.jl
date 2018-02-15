module Microbiome

export
    # Types
    DistanceMatrix,
    AbundanceTable,
    PCoA,

    # Functions
    ## abundance
    filterabund,
    relativeabundance,
    relativeabundance!,
    ## similarity
    getdm,
    getrowdm,
    pcoa,
    eigenvalue,
    principalcoord,
    optimalorder,
    optimalorder!,
    ## plotting
    hclustplot,
    annotationbar,
    ## utils
    panphlan_calcs,
    bysample,
    taxfilter,
    taxfilter!

using RecipesBase
using StatPlots
using StatsBase
using Distances
using Colors

import Clustering: Hclust, hclust
import DataFrames: DataFrame
import Base: getindex, setindex, length


include("utils.jl")
include("abundances.jl")
include("similarity.jl")
include("leafordering.jl")
include("plotting.jl")
include("biobakery_utils.jl")


end  # module Microbiome
