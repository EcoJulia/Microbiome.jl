module Microbiome

export
    # Types
    DistanceMatrix,
    AbundanceTable,
    PCoA,
        # re-exports
        Hclust,

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

        # re-exports
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


include("utils.jl")
include("abundances.jl")
include("similarity.jl")
include("plotting.jl")
include("biobakery_utils.jl")

end  # module Microbiome
