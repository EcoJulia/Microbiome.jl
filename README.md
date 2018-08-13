![Microbiome.jl logo](logo.png)

## For analysis of microbiome and microbial community data

**Latest Release:**

[![Latest Release](https://img.shields.io/github/release/BioJulia/Microbiome.jl.svg)](https://github.com/BioJulia/Microbiome.jl/releases/latest)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://biojulia.github.io/Microbiome.jl/stable)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/BioJulia/Microbiome.jl/blob/master/LICENSE)
![BioJulia maintainer: kescobo](https://img.shields.io/badge/BioJulia%20Maintainer-kescobo-blue.svg)
[![Build status](https://ci.appveyor.com/api/projects/status/wdpkdyafeadi5vx9?svg=true)](https://ci.appveyor.com/project/kescobo/microbiome-jl)


**Development Status**

[![Build Status](https://travis-ci.org/BioJulia/Microbiome.jl.svg?branch=master)](https://travis-ci.org/BioJulia/Microbiome.jl)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://biojulia.github.io/Microbiome.jl/latest)
[![Build status](https://ci.appveyor.com/api/projects/status/wdpkdyafeadi5vx9/branch/master?svg=true)](https://ci.appveyor.com/project/kescobo/microbiome-jl/branch/master)


## Description

Microbiome.jl is a package for manipulating and analyzing microbiome and
microbial community data. Many functions have been added to external packages
and are imported here.

## Installation

Install Microbiome from the Julia REPL:

```
pkg> add Microbiome
```
or

```julia
julia> using Pkg
julia> Pkg.add("Microbiome")
```

If you are interested in the cutting edge of the development, please check out
the master branch to try new features before release.

```
pkg> add Microbiome#master
```

or

```julia
julia> using Pkg
julia> Pkg.add("Microbiome#master")
```

## Companion Packages

You might also be interested in some functionality provided by
[MicrobiomePlots.jl](https://github.com/BioJulia/MicrobiomePlots)
and [BiobakeryUtils.jl](https://github.com/BioJulia/BiobakeryUtils).

Microbiome.jl uses [EcoBase](https://github.com/EcoJulia/EcoBase.jl) under the hood
for many of its types and methods.
