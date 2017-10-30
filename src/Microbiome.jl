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
    panphlan_calcs

using RecipesBase
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
