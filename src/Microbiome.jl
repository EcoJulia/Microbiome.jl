module Microbiome

export
    DistanceMatrix,
    AbundanceTable,

    # functions
    getdm

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

end  # module Microbiome
