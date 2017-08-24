module Microbiome

using BioSequences
using RecipesBase
using Distances
using DataFrames
using IterTools
using Distances

import MultivariateStats: classical_mds
import Clustering: Hclust, hclust
import Base: getindex, setindex, length

include("utils.jl")
include("abundances.jl")
include("similarity.jl")
include("plotting.jl")

end  # module Microbiome
