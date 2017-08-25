module Microbiome

export
    DistanceMatrix,
    AbundanceTable,

    # functions
    getdm

using BioSequences
using RecipesBase
using Distances
using DataFrames
using IterTools
using Distances

using MultivariateStats: classical_mds
using Clustering: Hclust, hclust
using Base: getindex, setindex, length

export DistanceMatrix,
    AbundanceTable,

    getdm


include("utils.jl")
include("abundances.jl")
include("similarity.jl")
include("plotting.jl")

end  # module Microbiome
