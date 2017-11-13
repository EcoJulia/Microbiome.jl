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
    relativeabundance!,
    hclustplot,
    panphlan_calcs,
    annotationbar

using RecipesBase
using StatPlots
using StatsBase
using Distances
using Colors

using DataFrames: DataFrame
using Clustering: Hclust, hclust
using Base: getindex, setindex, length


include("utils.jl")
include("abundances.jl")
include("similarity.jl")
include("plotting.jl")
include("biobakery_utils.jl")


end  # module Microbiome
