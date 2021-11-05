# News for Microbiome.jl

## v0.8.1

- Some changes for JOSS review ([#115](https://github.com/EcoJulia/Microbiome.jl/pull/115)),
  mostly in documentation
## v0.8.0

This release brings a host of changes,
including a manuscript for [submission to the Journal of Open Source Software](https://github.com/openjournals/joss-reviews/issues/3876).
Also, we're now part of [EcoJulia](http://github.com/EcoJulia)! (moved from [BioJulia](http://github.com/BioJulia)).

Notable changes include:

- Much better metadata handling in `CommunityProfile`s, (see esp issues [#63](https://github.com/EcoJulia/Microbiome.jl/issues/63) [#64](https://github.com/EcoJulia/Microbiome.jl/issues/64) [#65](https://github.com/EcoJulia/Microbiome.jl/issues/65) [#75](https://github.com/EcoJulia/Microbiome.jl/issues/75) [#76](https://github.com/EcoJulia/Microbiome.jl/issues/76))
- Add `commjoin` function to merge the contents of multiple `CommunityProfile`s ([#58](https://github.com/EcoJulia/Microbiome.jl/pull/58))
- Generic filtering of `CommunityProfile`s [#79](https://github.com/EcoJulia/Microbiome.jl/pull/79)

See [Release notes](https://github.com/EcoJulia/Microbiome.jl/releases/tag/v0.8.0) for full details.

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
