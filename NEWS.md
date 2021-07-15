# News for Microbiome.jl

## Unreleased

- Add `commjoin` function to merge the contents of multiple `CommunityProfile`s ([#58](https://github.com/BioJulia/Microbiome.jl/pull/58))

## v0.7.0 Major overhaul

Lots of breaking changes in this one, including

- Whole new community profile structure using `AxisArrays` back-end
- New `AbstractSample` type, containing embedded metadata
- New `Taxon` and `GeneFunction` feature types
- Dropped dependencies on `SpatialEcology.jl` and `DataFrames.jl`,
  but added `Tables.jl` interface and `EcoBase.jl` interfaces
- Brought back convenience `braycurtis` and `pcoa` functions relying on
  `Distances.jl` and `MultivariateStats.jl` respectively

## v0.5.0

Major Changes

- Dropped Microbiome.jl-defined types and methods for DistanceMatrix. Use Distances.jl instead.
  [Requires this SpatialEcology PR](https://github.com/EcoJulia/SpatialEcology.jl/pull/36)
- Dropped Microbiome.jl-defined types and methods for PCoA. Use MDS from MultivariateStats.jl instead. [Requires this MultivariateStats PR](https://github.com/JuliaStats/MultivariateStats.jl/pull/85)
- Dropped Hclust leaf ordering. Added to Clustering.jl instead ([see Clustering.jl PR](https://github.com/JuliaStats/Clustering.jl/pull/170))

## v0.4.1

Major Changes:

- Dropped support for julia-0.7

## v0.3.0

Major changes

- Dropped support for julia-0.6, moved on to julia-0.7
- Moved Plotting and BioBakery utilities out to separate packages to reduce dependencies
  - These two new packages will be released shortly
- I'm keeping DataFrames dependency for now, though I'd like to move to IterableTables soonish (though constructors might be more sense in EcoBase)

## Older versions

I didn't have NEWS.md until v0.3. I doubt anyone cares, but I kinda tried to
keep track in release notes.
