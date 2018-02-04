module Microbiome

export
    # Types
    DistanceMatrix,
    AbundanceTable,
    PCoA,

    # functions
    getdm,
    filterabund,
    getrowdm,
    pcoa,
    eigenvalue,
    principalcoord,
    relativeabundance,
    relativeabundance!,
    hclustplot,
    panphlan_calcs,
    annotationbar,
    optimalorder,
    optimalorder!,
    bysample

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
