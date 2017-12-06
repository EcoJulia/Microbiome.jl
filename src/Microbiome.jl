module Microbiome

export
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
    hclustplot,
    panphlan_calcs,
    annotationbar,

    # re-exports
    Hclust,
    DataFrame,
    hclust

using RecipesBase
using StatPlots
using StatsBase
using Distances
using Colors

import DataFrames: DataFrame
import Clustering: Hclust, hclust
import Base: getindex, setindex, length

using SpatialEcology

include("utils.jl")
include("abundances.jl")
include("similarity.jl")
include("plotting.jl")
include("biobakery_utils.jl")

end  # module Microbiome
