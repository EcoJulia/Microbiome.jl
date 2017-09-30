module Microbiome

export
    DistanceMatrix,
    AbundanceTable,
    PCoA,

    # functions
    getdm,
    filterabund,
    getrowdm,
    pcoa

using BioSequences
using RecipesBase
using Distances
using IterTools

using DataFrames: DataFrame
using Clustering: Hclust, hclust
using Base: getindex, setindex, length


include("utils.jl")
include("abundances.jl")
include("similarity.jl")
include("plotting.jl")
include("biobakery_utils.jl")


end  # module Microbiome
