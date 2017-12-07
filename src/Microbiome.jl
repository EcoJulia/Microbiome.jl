module Microbiome

export
    DistanceMatrix,
    abundancetable,
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

using Reexport
@reexport using SpatialEcology
import SpatialEcology.@forward_func


include("abundances.jl")
include("similarity.jl")
include("plotting.jl")
include("biobakery_utils.jl")

end  # module Microbiome
