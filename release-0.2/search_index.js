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
    "text": "Tables of abundances are based off ComMatrix types from SpatialEcology.jl, where columns are samples and rows are features (eg species). Sample and feature names are also stored, and there\'s a convenience function if you want to convert a DataFrame to a ComMatrix, assuming the first column contains feature names:using Microbiome\nusing DataFrames\n\ndf = DataFrame(species=[\"E. coli\", \"B. fragilis\", \"L. casei\"],\n                      sample1=[1, 4, 5],\n                      sample2=[3, 8, 0],\n                      sample3=[0, 3, 4]);\n\nabund = abundancetable(df)Forgive the clutter... ComMatricies name rows as species (which is true in this case, but need not be), and columns are \"sites\" rather than samples. That will be fixed eventually.samplenames(abund)\nfeaturenames(abund)\nsampletotals(abund) # this is column sums\nfeaturetotals(abund) # this is row sumsIf you want relative abundance, you can do relativeabundance(abund) or relativeabundance!(abund):relativeabundance!(abund)\n\nsampletotals(abund)You can also filter on the n most abundant features accross the dataset. This function automatically generates an n+1 row for other containing the remaining features. Note - these doesn\'t modify in-place, so you\'ve gotta reassign if you want to update:abund2 = filterabund(abund, 1)\n\nfeaturenames(abund2)"
},

{
    "location": "abundances.html#Plotting-1",
    "page": "Microbial Abundances",
    "title": "Plotting",
    "category": "section",
    "text": "Some convenience plotting types are available using RecipesBase and StatPlotsusing StatPlots\n\nsrand(1) # hide\n\nabund = abundancetable(\n           rand(100, 10),\n           [\"sample_$x\" for x in 1:10],\n           [\"feature_$x\" for x in 1:100]);\nrelativeabundance!(abund)\n\nabundanceplot(abund)\n\nsavefig(\"abundanceplot.png\"); nothing # hide(Image: )"
},

{
    "location": "distances.html#",
    "page": "Distances & Dissimilarity",
    "title": "Distances & Dissimilarity",
    "category": "page",
    "text": ""
},

{
    "location": "distances.html#Working-with-Distances-/-Dissimilarity-1",
    "page": "Distances & Dissimilarity",
    "title": "Working with Distances / Dissimilarity",
    "category": "section",
    "text": "Quite often, it\'s useful to boil stuff down to distances between samples. For this, I\'m using an interface with Distances.jl to generate a symetric DistanceMatrix, which also contains a vector for samples, and a field specifying which type of distance was used to calulate it. You can load one in manually, or generate it from an AbundanceTable.using Distances\nusing Microbiome\n\nabund = abundancetable([1  3  0;\n                        4  8  3;\n                        5  0  4]);\n\ndm = getdm(abund, BrayCurtis())I\'ve also implemented a method to do a principle coordinates analysis. If necessary, you can include correct_neg=true to use the correction method described in Lingoes (1971)p = pcoa(dm)\n\neigenvalue(p, 2)\nprincipalcoord(p, 1)\nvariance(p, [1,2])"
},

{
    "location": "distances.html#Plotting-1",
    "page": "Distances & Dissimilarity",
    "title": "Plotting",
    "category": "section",
    "text": "Some convenience plotting types are available using RecipesBase.using StatPlots\n\nsrand(1) # hide\nabund = abundancetable(\n    rand(100, 10),\n    [\"sample_$x\" for x in 1:10],\n    [\"feature_$x\" for x in 1:100]);\n\ndm = getdm(abund, BrayCurtis());\np = pcoa(dm, correct_neg=true);\n\nplot(p, title=\"Random PCoA\")\nsavefig(\"pcoplot.png\"); nothing # hide(Image: pcoa plot)"
},

{
    "location": "distances.html#Optimal-Leaf-Ordering-1",
    "page": "Distances & Dissimilarity",
    "title": "Optimal Leaf Ordering",
    "category": "section",
    "text": "I\'ve also provided a plotting recipe for making treeplots for Hclust objects from the Clustering.jl package:using Clustering\n\ndm = [0. .1 .2\n      .1 0. .15\n      .2 .15 0.];\n\nh = hclust(dm, :single);\nh.labels = [\"a\", \"b\", \"c\"];\n\nhclustplot(h)\nsavefig(\"hclustplot1.png\"); nothing # hide(Image: hclust plot 1)Note that even though this is a valid tree, the leaf a is closer to leaf c, despite the fact that c is more similar to b than to a. This can be fixed with a method derived from the paper:Bar-Joseph et. al. \"Fast optimal leaf ordering for hierarchical clustering.\" _Bioinformatics_. (2001)optimalorder!(h, dm)\nhclustplot(h)\n\nsavefig(\"hclustplot2.png\"); nothing # hide(Image: hclust plot 1)"
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
    "text": "The BioJulia organisation has a set of contribution guidelines which apply to all BioJulia projects. These guidelines are available here and it is recommended all new contributors read these guidelines before opening a pull request."
},

{
    "location": "contributing.html#Making-a-contribution-1",
    "page": "Contributing",
    "title": "Making a contribution",
    "category": "section",
    "text": "If you\'re interested in adding functionality to Microbiome.jl, please feel free to open an issue or a pull request (PR) against the master branch. If you\'re not yet ready for that, you can also ask questions/start a discussion in the Bio.jl gitter channel. Work-in-progress PRs are fine, as discussion about approach and code review can happen in the PR.Before merging, any new code should be unit tested and have docs for newly exported functions, but if you don\'t know how to do this, don\'t worry, we can help!"
},

]}
