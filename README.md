![Microbiome.jl logo](logo.png)

[![status](https://joss.theoj.org/papers/450fa18f47932c5fd3b837edeac91440/status.svg)](https://joss.theoj.org/papers/450fa18f47932c5fd3b837edeac91440) ![EcoJulia maintainer: kescobo](https://img.shields.io/badge/EcoJulia%20Maintainer-kescobo-blue.svg)
## Latest Release 

[![Latest Release](https://img.shields.io/github/release/EcoJulia/Microbiome.jl.svg)](https://github.com/EcoJulia/Microbiome.jl/releases/latest)

[![Docs stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://docs.ecojulia.org/Microbiome.jl/stable/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/EcoJulia/Microbiome.jl/blob/master/LICENSE)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

## Development Status

[![Docs dev](https://img.shields.io/badge/docs-latest-blue.svg)](https://docs.ecojulia.org/Microbiome.jl/latest/)
[![CI](https://github.com/EcoJulia/Microbiome.jl/workflows/CI/badge.svg)](https://github.com/EcoJulia/Microbiome.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/EcoJulia/Microbiome.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/EcoJulia/Microbiome.jl)

## Description

Microbiome.jl is a package for manipulating and analyzing
microbiome and microbial community data.

## Installation

To use the latest version of `Microbiome.jl`,
you must be on `julia` v1.6 or greater.

Install `Microbiome.jl` from the Julia REPL:

```julia
julia> using Pkg

julia> pkg"add Microbiome"
```

If you are interested in the cutting edge of the development, please check out
the master branch to try new features before release.

```julia
julia> pkg"add Microbiome#main"
```

## Companion Packages

- You might also be interested in some functionality provided by
  [BiobakeryUtils.jl](https://github.com/EcoJulia/BiobakeryUtils.jl).
- Microbiome.jl uses [EcoBase](https://github.com/EcoJulia/EcoBase.jl) under the hood
  for many of its types and methods.
- Microbiome.jl has a `Tables.jl` interface, so you can convert `CommunityProfiles`
  to any `Tables.jl`-compatible sink.
  - `CommunityProfile` is not yet a sink itself
