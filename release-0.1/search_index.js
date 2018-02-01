var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Microbiome.jl-smallFor-analysis-of-microbiome-and-microbial-community-data/small-1",
    "page": "Home",
    "title": "Microbiome.jl <small>For analysis of microbiome and microbial community data</small>",
    "category": "section",
    "text": "(Image: Latest Release) (Image: Microbiome) (Image: License) (Image: ) (Image: ) (Image: BioJulia maintainer: kescobo)Development builds: (Image: Build Status)"
},

{
    "location": "index.html#Description-1",
    "page": "Home",
    "title": "Description",
    "category": "section",
    "text": "Microbiome.jl is a package for manipulating and analyzing microbiome and microbial community data. Many functions have been added to external packages and are imported here."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "Install Microbiome from the Julia REPL:julia> Pkg.add(\"Microbiome\")If you are interested in the cutting edge of the development, please check out the master branch to try new features before release.julia> Pkg.checkout(\"Microbiome\")"
},

{
    "location": "abundances.html#",
    "page": "Microbial Abundances",
    "title": "Microbial Abundances",
    "category": "page",
    "text": ""
},

{
    "location": "abundances.html#Working-with-microbial-abundances-1",
    "page": "Microbial Abundances",
    "title": "Working with microbial abundances",
    "category": "section",
    "text": "The AbundanceTable type is treated like a 2D array where columns are samples and rows are features (eg species). Sample and feature names are also stored, and there's a convenience function if you want to convert a DataFrame to an AbundanceTable, assuming the first column contains feature names:julia> using Microbiome\njulia> using DataFrames\n\njulia> df = DataFrame(species=[\"E. coli\", \"B. fragilis\", \"L. casei\"],\n                      sample1=[1, 4, 5],\n                      sample2=[3, 8, 0],\n                      sample3=[0, 3, 4])\n3×4 DataFrames.DataFrame\n│ Row │ species     │ sample1 │ sample2 │ sample3 │\n├─────┼─────────────┼─────────┼─────────┼─────────┤\n│ 1   │ E. coli     │ 1       │ 3       │ 0       │\n│ 2   │ B. fragilis │ 4       │ 8       │ 3       │\n│ 3   │ L. casei    │ 5       │ 0       │ 4       │\n\njulia> abund = AbundanceTable(df)\n3×3 Microbiome.AbundanceTable{Int64}:\n 1  3  0\n 4  8  3\n 5  0  4If you want relative abundance, you can do relativeabundance(abund) or relativeabundance!(abund):julia> abund = relativeabundance(abund)\n3×3 Microbiome.AbundanceTable{Float64}:\n 0.1  0.272727  0.0\n 0.4  0.727273  0.428571\n 0.5  0.0       0.571429\n ```\n\nYou can also filter on the `n` most abundant features accross the dataset. This\nfunction automatically generates an `n+1` row for `other` containing the\nremaining features. Note - these doesn't modify in-place, so you've gotta\nreassign if you want to update:\njulia julia> abund2 = filterabund(abund, 1) 2×3 Microbiome.AbundanceTable{Float64}:  0.4  0.727273  0.428571  0.6  0.272727  0.571429julia> abund2.features 2-element Array{String,1}:  \"B. fragilis\"  \"other\"\n## Plotting\n\nSome convenience plotting types are available using [`RecipesBase`][1] and\n[StatPlots][2]\n\n[1]: https://github.com/juliaplots/recipesbase.jl\n[2]: https://github.com/juliaplots/StatPlots.jl\njulia using StatPlots # TODO: add actual exampleabund = AbundanceTable(     rand(100, 10), [\"sample_x for x in 110     feature_x\" for x in 1:100])abund = relativeabundance(abund) plot(abund, title=\"Random abundance\") ```"
},

{
    "location": "similarities.html#",
    "page": "Distances & Dissimilarity",
    "title": "Distances & Dissimilarity",
    "category": "page",
    "text": ""
},

{
    "location": "similarities.html#Working-with-Distances-/-Dissimilarity-1",
    "page": "Distances & Dissimilarity",
    "title": "Working with Distances / Dissimilarity",
    "category": "section",
    "text": "Quite often, it's useful to boil stuff down to distances between samples. For this, I'm using an interface with Distances.jl to generate a symetric DistanceMatrix, which also contains a vector for samples, and a field specifying which type of distance was used to calulate it. You can load one in manually, or generate it from an AbundanceTable.julia> using Distances\njulia> using Microbiome\n\njulia> abund = AbundanceTable([1  3  0;\n                               4  8  3;\n                               5  0  4])\n3×3 Microbiome.AbundanceTable{Int64}:\n 1  3  0\n 4  8  3\n 5  0  4\n\njulia> dm = getdm(abund, BrayCurtis())\n3×3 Microbiome.DistanceMatrix{Float64}:\n 0.0       0.52381   0.176471\n 0.52381   0.0       0.666667\n 0.176471  0.666667  0.0I've also implemented a method to do a principle coordinates analysis. If necessary, you can include correct_neg=true to use the correction method described in [Lingoes (1971)][2]julia> p = pcoa(dm)\n\n3×2 Microbiome.PCoA{Float64}:\n -0.251198   0.776895\n  0.79841   -0.170903\n -0.547212  -0.605992\n\njulia> eigenvalue(p, 2)\n0.0050620487880063784\n\njulia> principalcoord(p, 1)\n3-element Array{Float64,1}:\n -0.251198\n  0.79841\n -0.547212\n\njulia> p.variance_explained\n2-element Array{Float64,1}:\n 0.979751\n 0.0202492"
},

{
    "location": "similarities.html#Plotting-1",
    "page": "Distances & Dissimilarity",
    "title": "Plotting",
    "category": "section",
    "text": "Some convenience plotting types are available using [RecipesBase][1].[1]: https://github.com/juliaplots/recipesbase.jlusing StatPlots # TODO: add actual example\n\nabund = AbundanceTable(\n    rand(100, 10), [\"sample_$x\" for x in 1:10],\n    [\"feature_$x\" for x in 1:100])\n\ndm = getdm(abund, BrayCurtis())\np = pcoa(dm, correct_neg=true)\n\nplot(p, title=\"Random PCoA\")"
},

{
    "location": "contributing.html#",
    "page": "Contributing",
    "title": "Contributing",
    "category": "page",
    "text": ""
},

{
    "location": "contributing.html#Contributing-1",
    "page": "Contributing",
    "title": "Contributing",
    "category": "section",
    "text": "The BioJulia organisation has a set of contribution guidelines which apply to all BioJulia projects. These guidelines are available [here][1] and it is recommended all new contributors read these guidelines before opening a pull request.[1]: http://biojulia.github.io/Contributing/latest"
},

{
    "location": "contributing.html#Making-a-contribution-1",
    "page": "Contributing",
    "title": "Making a contribution",
    "category": "section",
    "text": "If you're interested in adding functionality to Microbiome.jl, please feel free to open an issue or a pull request (PR) against the master branch. If you're not yet ready for that, you can also ask questions/start a discussion in the [Bio.jl][2] gitter channel. [2]: https://gitter.im/BioJulia/Bio.jlWork-in-progress PRs are fine, as discussion about approach and code review can happen in the PR.Before merging, any new code should be unit tested and have docs for newly exported functions, but if you don't know how to do this, don't worry, we can help!"
},

]}
