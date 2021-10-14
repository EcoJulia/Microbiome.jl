---
title: 'Microbiome.jl and BiobakeryUtils.jl - Julia packages for working with microbial community data'
tags:
  - julia
  - biology
  - microbiome
  - microbial communities
authors:
  - name: Kevin S. Bonham^[Corresponding author]
    affiliation: 1
  - name: Annelle Abatoni Kayisire
    affiliation: 1
  - name: Anika Luo
    affiliation: 1
  - name: Vanja Klepac-Ceraj^[Corresponding author]
    affiliation: 1
affiliations:
 - name: Department of Biological Sciences, Wellesley College
   index: 1
date: 30 October 2021
bibliography: paper.bib
---

# Summary

`Microbiome.jl` is a julia package to facilitate analysis of microbial community data.
`BiobakeryUtils.jl` is built on top of `Microbiome.jl`,
and provides utilities for working with a suite of command line tools
(the bioBakery) that are widely used for converting raw metagenomic sequencing data
into tables of taxon and gene function counts.
Together, these packages provide an effective way to link microbial community data
with the power of julia's numerical, statistical, and plotting libraries.

# Statement of need

Complex microbial communities exist everywhere, including in and on the human body,
and have profound effects on the environment and human health [@LloydPrice2017].
Common methods for analyzing microbial communities (eg 16S amplicon or metagenomic sequencing)
generate a large quantity of numerical data (eg count or relative abundance data)
as well as metadata associated with biological samples (locations, human subject data)
and microbial features (taxa, gene functions) [@Mallick2017ExperimentalDA].

The julia programming language [@Bezanson2017-ud] is gaining increasing prominence in biological research
due to its speed and flexibility [@roesch2021julia],
and has a growing ecosystem of packages for working with biological and ecological data,
as well as libraries for Bayesian statistical analysis [@ge2018t],
scientific machine learning [@rackauckas2017differentialequations],
and plotting [@DanischKrumbiegel2021].
Julia's type system makes it incredibly easy for packages to interoperate,
making `Microbiome.jl` and `BiobakeryUtils.jl` and effective bridge between
microbial community data and julia's package ecosystem,
while remaining agnostic to downstream analysis.

# Functionality

# Limitations and future work
# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures


# Acknowledgements

<!-- TODO -->

# References