### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 38e1904d-b7bb-4868-af51-38688c8f65e6
using Downloads

# ╔═╡ cfc9a88f-4b7e-48e6-800b-830e7d368597
using CSV, DataFrames, Microbiome

# ╔═╡ d5390ec3-c8c0-4025-92ac-42c762ec361b
using Distances, MultivariateStats, Plots

# ╔═╡ bd9652fc-f262-11ec-3cf6-71de5b26f0f0
md"""
# 🦠 Microbiome.jl tutorial

In this tutorial,
we'll learn how to interact with microbial community data
using [`Microbiome.jl`](https://github.com/EcoJulia/Microbiome.jl).
If you have any questions or comments,
please start a discussion
or open an issue
on github!

Let's get started! 👇️
"""

# ╔═╡ 8d01fe33-6610-44b1-800a-5fba9224359c
md"""
## ⬇️ Getting some example data

The first thing we'll need is some data!
If you have your own data, you should be able to substitute it below.
If you don't, we can use some example data
from the human microbiome project
(here downloaded from the [`MetaPhlan` tutorial](https://github.com/biobakery/biobakery/wiki/metaphlan3#run-multiple-samples)).

These files are subsampled from real human biosamples
from different body sites:

- SRS014459-Stool_profile.txt (💩)
- SRS014464-Anterior_nares_profile.txt (back of the 👃)
- SRS014470-Tongue_dorsum_profile.txt (top of the 👅)
- SRS014472-Buccal_mucosa_profile.txt (inside of cheek)
- SRS014476-Supragingival_plaque_profile.txt (🦷 plaque near the gums)
- SRS014494-Posterior_fornix_profile.txt (vagina near the cervix)
"""

# ╔═╡ 1be1644b-7270-4a8d-a208-8f000df47c5d
files = [
	"SRS014459-Stool_profile.txt",
	"SRS014464-Anterior_nares_profile.txt",
	"SRS014470-Tongue_dorsum_profile.txt",
	"SRS014472-Buccal_mucosa_profile.txt",
	"SRS014476-Supragingival_plaque_profile.txt",
	"SRS014494-Posterior_fornix_profile.txt"
]

# ╔═╡ 4cf6bb90-684b-431e-ae72-fdf541550d17
baseurl = "https://github.com/biobakery/biobakery/raw/master/test_suite/biobakery_tests/data/metaphlan3/output/";

# ╔═╡ c78b8938-2f4e-4b21-afa0-01c6aa03b689
for f in files
	localpath = joinpath(".", f)
	isfile(localpath) && continue # skip file if it's already downloaded
	Downloads.download(joinpath(baseurl, f), localpath)
end

# ╔═╡ 64624649-8a70-44fb-bbd5-5c3466ce484a
@assert all(isfile, files)

# ╔═╡ ce747551-460b-4542-b07c-7adc3fdd08cd
md"""
### 📊 Data format

The data we've downloade here is in the form of tab-separated files,
one per biosample,
where rows represent individual features - in this case, indivudual taxa.

Let's take a look:
"""

# ╔═╡ 8a45213f-2429-4035-a782-81e5bf8fba17
first(eachline("SRS014459-Stool_profile.txt"), 6)

# ╔═╡ 37940aa2-2876-4d3a-bed4-306aa42d5f42
md"""
The first few rows are headers showing how the file was generated,
the actual data starts on line 4.
The first and second columns contain the taxonomy of the feature
(in human readable form and in NCBI's ID)
and the third column
is the relative abundance of that taxon.
The final column contains information about any additional species
that may be contained in a given identified row,
but we can ignore that for now.

Let's load it into a DataFrame to start with:
"""

# ╔═╡ 370e4f43-3835-45d1-b614-5732c66842f5
stool = CSV.read("SRS014459-Stool_profile.txt", DataFrame;
	header = 4, # get the column headers from the 4th line
	skipto = 5, # start parsing on the 5th line
	delim = '\t', # use tabs 
)

# ╔═╡ 88564be4-7943-4577-adee-b5e96351912b
md"""
## The `CommunityProfile` type

The main contribution of `Microbiome.jl` is the `CommunityProfile` type,
which consists of 3 parts:

1. A sparse matrix containing the relative abundances of features.
2. A list of features. These index the rows of the matrix
3. A list of biosamples with optionally attached metadata
   (see `MicrobiomeSample` below). These index the columns of the matrix

We'll contruct it mannually first,
and then later we'll see how we can use the `Tables.jl` integration
to our advantage.
But first, let's load in all of our samples:
"""

# ╔═╡ 4a5cfe7c-109f-4bb0-ba17-6a4aba9dc16a
comm_df = let dfs = DataFrame[]
	for f in files
		df = CSV.read(f, DataFrame;
			header = 4, # get the column headers from the 4th line
			skipto = 5, # start parsing on the 5th line
			delim = '\t', # use tabs 
		)
	
		samplename=splitext(f)[1] # get sample name from file
		
		rename!(df, "#clade_name" => "taxon", "relative_abundance" => samplename)
		push!(dfs, select(df, [1,3]))
	end
	outerjoin(dfs...; on="taxon")
end

# ╔═╡ 287fd231-4ff9-4c2e-a55d-1786a7e41659
md"""
### 🏁 Relative abundance matrix

You can use something from the `SparseArrays` standard library,
but any old matrix will do, and `Microbiome.jl` will convert it for you.

Notice that, in the DataFrame above,
a number of taxa are missing from a number of samples.
This is expected,
but "missing" in the contex of microbial community data just means
that the relative abundance is 0.
So we'll convert those to 0 to generate our matrix:
"""

# ╔═╡ 851b7260-1701-4e92-b2ec-cc9ba37d2e46
mat = coalesce.( # this function returns either a non-missing value or 0
	Matrix(comm_df[:, 2:end]), # our numerical values are in columns 2:end
	0.0
)

# ╔═╡ cdabefa8-77c6-4a56-bef4-c5591330ea59
md"""
### Microbiome Features (taxa, gene functions, metabolites)

The rows in a `CommunityProfile` are instances of `AbstractFeature`s.
`AbstractFeature`s need have only a `name` field,
but there are some concrete types provided by `Microbiome.jl`
that have some additional properties.
In this case, because we're looking at taxonomic profiles,
we'll use the `Taxon` type, which also contains a `rank` field
that tells you the taxonomic rank (eg Kingdom, Phylum, Genus, etc).

You may or may not have noticed that our table in many cases has multiple rows
delineating different taxonomic ranks.
"""

# ╔═╡ 425cb918-4a8b-46e2-810b-17d32c7081f5
first(comm_df.taxon, 5)

# ╔═╡ 2908dc16-74ee-4c87-acf8-681b3f2854f2
md"""
Here, the order "Lactobacillales" is also represented
in the class "Bacilli" and the phylum "Firmicutes" etc.

To make our `Taxon` type, we could manually extract
the rank and name of the last element...
"""

# ╔═╡ 239406a6-0424-4107-a809-2fccb1dfb8b1
lacto_ranks = split(comm_df.taxon[4], '|')

# ╔═╡ 1fb80604-2385-49d6-b07c-b6b9f467c616
# get the components of the last element
lacto_rank, lacto_name = split(last(lacto_ranks), "__")

# ╔═╡ 5456d0d6-2b69-4a91-9ef2-0ed55361c39f
Taxon(lacto_name, :order)

# ╔═╡ bde7469f-4734-4882-bd3b-2971468dfa3e
md"""
... or we can use the `taxon()` function,
which understands this commonly used format and will build it for us:
"""

# ╔═╡ a42a53d3-8046-4c2a-bbef-4cd4d246b141
taxon(last(lacto_ranks))

# ╔═╡ 9812b45d-4973-414e-a71c-22dd81367f88
md"""
So, to get a vector of `Taxon`s, we can do the following:
"""

# ╔═╡ 064c6b27-f16d-4fbd-9ec6-67e091ecc9d2
taxa = map(comm_df.taxon) do t
	ranks = split(t, '|')
	taxon(last(ranks))
end

# ╔═╡ 4e282e00-0749-4c0a-b44d-ef3cdfbdb10f
md"""
### The `MicrobiomeSample` type

Biosamples in `Microbiome.jl` index the columns of our abundance matrix,
and can also store optional metadata.
Let's get the sample names from our file names,
and then add the body site as metadata.

We'll use the first sample as an example:
"""

# ╔═╡ 83ebfa19-b6d6-40db-9cae-bb5e840b4e46
stool_file = first(files)

# ╔═╡ 36591c12-f038-4ea7-af94-1e60b7aad0ba
(samplename, site, ext) = split(stool_file, r"(-|_profile)"; limit = 3)

# ╔═╡ 0af2ae11-cc91-44a2-8fb9-cef16d15fbee
stool_samp = MicrobiomeSample(samplename)

# ╔═╡ df70cfe0-2859-4b22-ba01-ed12db675e27
md"""
Sample metadata is stored in a `Dictionary` from the `Dictionaries.jl` package.
We can add metadata using `set!` or `insert!`
(the latter will throw an error if the metadata already exists):
"""

# ╔═╡ 4de0df7b-457d-4a38-a47f-83e536c47473
set!(stool_samp, :body_site, site)

# ╔═╡ 8715a036-c3b3-49bd-af31-d349c0d63bc3
md"""
So let's do this for all of our samples:
"""

# ╔═╡ b9d102ac-0aab-45f6-9f1d-f23ba590bcc6
mbiom_samples = map(files) do f
	(name, site, ext) = split(f, r"(-|_profile)"; limit = 3)
	s = MicrobiomeSample(name)
	insert!(s, :body_site, site)
	s
end

# ╔═╡ 5e77d309-ca1b-47b2-8c71-a0cb2c96758a
md"""
### Constructing a `CommunityProfile`

Now have all the elements, and we just need to put them together:
"""

# ╔═╡ 96ca485b-db0a-4517-b15b-6cb3232c86e2
comm = CommunityProfile(mat, taxa, mbiom_samples)

# ╔═╡ f5decd43-03ca-40ce-89a0-2cb74fa0fa70
md"""
## Working with `CommunityProfile`s

Now that we have our `CommunityProfile`, there are a number of ways to work with it.
First, let's take a look at the various accessor functions.

- `features()` - returns a vector of the features (in this case `Taxon`s)
- `samples()` - returns a vector of the samples, including the attached metadata
- `featurenames()`/`samplenames()` - like `features()` and `samples()`, but returns just the string representation.
- `abundances()` - returns the sparse matrix with relative abundances`
- `ranks()` - if features are `Taxon`s, returns a vector of the ranks of those taxa (if they have them)
"""

# ╔═╡ 5c85e0d5-d0aa-4e32-8ecd-e2688fc51201
features(comm)

# ╔═╡ 67bc83af-a427-4380-81fb-dbdb14cd2439
samples(comm)

# ╔═╡ 58a2aa85-8cce-4986-bfa2-5284afda40e6
featurenames(comm)

# ╔═╡ de790428-9248-4ec3-8bb6-b8fd187483ec
ranks(comm)

# ╔═╡ 95011cf6-e81b-4a48-9fdc-388df2c9ed75
md"""
### ℹ️ Indexing

There are a number of useful ways to index into a `CommunityProfile`.
The first and most straightforward is using numbers -
keep in mind that indexing in julia starts at 1,
that rows are the first dimension, and columns are the second dimension.

For example, the following will return the relative abundance of the 3rd row and 6th column:
"""

# ╔═╡ c98a4d00-e94d-47aa-aa1b-d91c8422ddee
comm[4, 6]

# ╔═╡ 44e657e2-01bf-4fea-a2c3-c9eb900a6214
md"""
One can also get slices using standard julia syntax.
Note that when taking slices, you get back another `CommmunityProfile`.
"""

# ╔═╡ e1f01522-8544-47ca-94ba-0deb9acaba15
comm[1:3, 2:4]

# ╔═╡ d3c9d708-f556-4ef6-a34b-c7e461351af4
md"""
One huge benefit of a `CommunityProfile`
is that we can also index using strings and regular expressions
that match on the `name` field of `AbstractFeature`s (rows)
and `MicrobiomeSample`s:
"""

# ╔═╡ 9f2b2856-4256-48dc-8f02-a45dbc898798
comm[r"Baci", :]

# ╔═╡ ad1cec80-6f26-4541-a6af-2363257f609d
md"""
### 🌳 Ecological metrics

There are a couple of built-in functions for calculating
things like alpha diversity of samples.
For example, to caclutate the Gini-Simpson or Shannon indecies,
use `ginisimpson()` or `shannon()` respectively, both of which
return a 1 x N matrix with length equal to the number of samples you have:
"""

# ╔═╡ 8529552b-bb06-49bd-937a-0bf21059b144
shannon(comm)

# ╔═╡ 15a5f712-415b-41c4-a3b0-d236075a7aba
md"""
Alternatively, you can use the mutating versions of these functions
(`ginisimpson!()` and `shannon!()`) to automatically add the result
to the metadata for each sample.
"""

# ╔═╡ c641f018-5e61-407c-a651-7c7f7b114449
ginisimpson!(comm);

# ╔═╡ 24aa8793-5470-4122-acd7-92a29312689b
shannon!(comm);

# ╔═╡ f322ff0a-49e4-41a2-96d7-e73d5c1f29bd
samples(comm)

# ╔═╡ 351a621b-c782-45dd-b141-d248eb07e8e8
md"""
We can also get information about the presence and prevalence of indivudual features.
The `present()` function takes a `CommunityProfile`
and an optional minimum abundance to be considered present (default > 0)
and returns a boolean matrix:
"""

# ╔═╡ 8b6bc20b-6cc6-4450-a8c3-9b58fdbda6d8
present(comm, 1.0) # anything over 1%

# ╔═╡ 1e48e79b-7242-473b-a730-0ea7c5313c0f
md"""
The `prevalence()` function calculates the fraction of samples
that a given feature is found, also with an optional minimum abundance value:
"""

# ╔═╡ 0b80a1b8-b6c6-4663-987d-b4613536056b
prevalence(comm, 1.0)

# ╔═╡ 7fdf6e5c-471e-4170-9a7c-2c027ec02100
md"""
## Working with sample metadata

One of the major benefits of using `Microbiome.jl` is coupling metadata
to samples, enabling simple downstream analysis.

We've already seen how you can add individual metadata values to `MicrobiomeSample`s.
But you can also do this when the samples are embedded in ComminityProfiles.
"""

# ╔═╡ cd06d809-0a85-485e-9ec7-7698cc748344
set!(comm, "SRS014459", :new_metadatum, "some value");

# ╔═╡ 3a064eac-2e45-4565-a35b-71e3f755c5c0
samples(comm)[1]

# ╔═╡ e57c4876-f4c2-40d4-befc-1f6e82ae1496
md"""
Or, add lots of metadata to your community profile
using any `Tables.jl` compatible table.
All that's required is a column called "sample" to match metadata rows.
"""

# ╔═╡ e5cbfd13-7224-453d-b108-a4e54259540d
mtdt = DataFrame(
	sample = samplenames(comm),
	sex    = rand(("male", "female"), 6),
	age    = [10, 11, 12, 13, 14, 15]
)

# ╔═╡ bc0f81b3-3750-426e-8113-712f7490236e
set!(comm, mtdt)

# ╔═╡ f151ec8d-08c4-4e2c-bc2e-555b999d771a
md"""
And you can get all of the metadata out as a vector of `NamedTuple`s
(which can be converted into `DataFrame`s or other Table types).
"""

# ╔═╡ 71be47ae-1682-4068-a4d7-e36291682efa
DataFrame(metadata(comm))

# ╔═╡ 8d64359d-3ab4-4da0-af92-afc4d25389bf
md"""
Notice that the `"new_metadatum"` column gets filled with `missing`
for the samples that don't have an entry.

You can also get out a vector of individual metadata values using `get()`:
"""

# ╔═╡ e0675c3d-7e3d-4866-adc0-ade649d3955c
get(comm, :body_site)

# ╔═╡ 42ff9700-8d33-41b0-915a-acb6202bc227
md"""
## Integration with other Julia packages

Now let's see how we can use the features of `CommunityProfile`s
to do an analysis of microbiome data.
Because the underlying data is just a matrix,
We can integrate out microbial community data with the wide array
of julia packages for doing frequentist or Bayesian statistics,
machine learning, and plotting.

As a brief example, let's fit a principal coordinates analysis (PCoA)
model and plot it using some of our metadata.

NOTE: usually you would first filter on a single taxonmic rank,
but since the example data is so sparse, it's really boring when we do that.
"""

# ╔═╡ 82d5d99b-e2b6-46a0-806e-514a019ccb4d
#eg
filter(t-> taxrank(t) == :species, comm)

# ╔═╡ 72a0d2de-aa12-492f-a4f2-22b5eb8e5076
# the `pairwise()` function is from `Distances.jl`
distance_matrix = pairwise(BrayCurtis(), abundances(comm); dims=2)

# ╔═╡ bee631cc-5aa4-45e5-9499-38e60bad9418
# from `MultivariateStats.jl`
pco = fit(MDS, distance_matrix; distances=true)

# ╔═╡ 184c8e04-687b-4ea7-90ad-40fb12cd5b13
# each column is a principal coordinates axis
pco_axes = projection(pco)

# ╔═╡ e59580c0-8b15-417a-951c-0e7778a949d9
scatter(pco_axes[:, 1], pco_axes[:, 2]; 
	xlabel = "PCo Axis 1",
	ylabel = "PCo Axis 2",
	legend = false
)

# ╔═╡ 0eabcaa7-4c33-4463-a3bb-c4da1f37b36a
scatter(pco_axes[:, 1], pco_axes[:, 2]; 
	xlabel = "PCo Axis 1",
	ylabel = "PCo Axis 2",
	title  = "Colored by body site",
	groups = get(comm, :body_site),
	legend = :bottomright
)

# ╔═╡ 9eb41aa2-f534-488a-8f1b-90842353393d
scatter(pco_axes[:, 1], pco_axes[:, 2]; 
	xlabel = "PCo Axis 1",
	ylabel = "PCo Axis 2",
	title  = "Colored by Shannon diversity", 
	zcolor = get(comm, :shannon),
	legend = :bottomright
)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Distances = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
Downloads = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
Microbiome = "3bd8f0ae-a0f2-5238-a5af-e1b399a4940c"
MultivariateStats = "6f286f6a-111f-5878-ab1e-185364afe411"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"

[compat]
CSV = "~0.10.4"
DataFrames = "~1.3.4"
Distances = "~0.10.7"
Microbiome = "~0.9.1"
MultivariateStats = "~0.9.1"
Plots = "~1.31.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra", "Logging"]
git-tree-sha1 = "91ca22c4b8437da89b030f08d71db55a379ce958"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.5.3"

[[deps.Arpack_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "5ba6c757e8feccf03a1554dfaf3e26b3cfc7fd5e"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.1+1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "873fb188a4b9d76549b81465b1f75c82aaf59238"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.4"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9489214b993cd42d17f44c36e359bf6a7c919abf"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "1e315e3f4b0b7ce40feded39c73049692126cf53"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.3"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "1fd869cc3875b57347f7027521f561cf46d1fcd8"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.19.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "9be8be1d8a6f44b96482c8af52238ea7987da3e3"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.45.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "daa21eb85147f72e41f6352a57fccea377e310a9"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.4"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Dictionaries]]
deps = ["Indexing", "Random"]
git-tree-sha1 = "7669d53b75e9f9e2fa32d5215cb2af348b2c13e2"
uuid = "85a47980-9c8c-11e8-2b9f-f7ca1fa99fb4"
version = "0.3.21"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.EcoBase]]
deps = ["RecipesBase"]
git-tree-sha1 = "a4d5b263972e820e780effc2084f92399ba44ee3"
uuid = "a58aae7d-b440-5a11-b283-399458f99aac"
version = "0.1.6"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "129b104185df66e408edd6625d480b7f9e9823a0"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.18"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "51d2dfe8e590fbd74e7a842cf6d13d8a2f45dc01"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.6+0"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "c98aea696662d09e215ef7cda5296024a9646c75"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.64.4"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "3a233eeeb2ca45842fe100e0413936834215abf5"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.64.4+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "83ea630384a13fc4f002b77690bc0afeb4255ac9"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.2"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Indexing]]
git-tree-sha1 = "ce1566720fd6b19ff3411404d4b977acd4814f9f"
uuid = "313cdc1a-70c2-5d6a-ae34-0150d3930a38"
version = "1.1.1"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "61feba885fac3a407465726d0c330b3055df897f"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.2"

[[deps.InlineTest]]
deps = ["Test"]
git-tree-sha1 = "daf0743879904f0ad645ca6594e1479685f158a2"
uuid = "bd334432-b1e7-49c7-a2dc-dd9149e4ebd6"
version = "0.2.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "b3364212fb5d870f724876ffcd34dd8ec6d98918"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.7"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "46a39b9c58749eefb5f2dc1178cb8fab5332b1ab"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.15"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "09e4b894ce6a976c354a69041a04748180d43637"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.15"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Microbiome]]
deps = ["Dictionaries", "Distances", "EcoBase", "MultivariateStats", "ReTest", "SparseArrays", "Statistics", "Tables"]
git-tree-sha1 = "e3393084d259a6c130e85461de1bad059967bc61"
uuid = "3bd8f0ae-a0f2-5238-a5af-e1b399a4940c"
version = "0.9.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.MultivariateStats]]
deps = ["Arpack", "LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI", "StatsBase"]
git-tree-sha1 = "7008a3412d823e29d370ddc77411d593bd8a3d03"
uuid = "6f286f6a-111f-5878-ab1e-185364afe411"
version = "0.9.1"

[[deps.NaNMath]]
git-tree-sha1 = "737a5957f387b17e74d4ad2f440eb330b39a62c5"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9a36165cf84cff35851809a40a928e1103702013"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.16+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "0044b23da09b5608b4ecacb4e5e6c6332f833a7e"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.2"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "8162b2f8547bc23876edd0c5181b27702ae58dce"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.0.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "9888e59493658e476d3073f1ce24348bdc086660"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.0"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "93e82cebd5b25eb33068570e3f63a86be16955be"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.31.1"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "c6c0f690d0cc7caddb74cef7aa847b824a16b256"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.ReTest]]
deps = ["Distributed", "InlineTest", "Printf", "Random", "Sockets", "Test"]
git-tree-sha1 = "dd8f6587c0abac44bcec2e42f0aeddb73550c0ec"
uuid = "e0db7c4e-2690-44b9-bad6-7687da720f89"
version = "0.3.2"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "dc1e451e15d90347a7decc4221842a022b011714"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.2"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "db8481cf5d6278a121184809e9eb1628943c7704"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.13"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "a9e798cae4867e3a41cae2dd9eb60c047f1212db"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.6"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "9f8a5dc5944dc7fbbe6eb4180660935653b0a9d9"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.0"

[[deps.StaticArraysCore]]
git-tree-sha1 = "66fe9eb253f910fe8cf161953880cfdaef01cdf0"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.0.1"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "8d7530a38dbd2c397be7ddd01a424e4f411dcc41"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.2"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "642f08bf9ff9e39ccc7b710b2eb9a24971b52b1a"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.17"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "ec47fb6069c57f1cee2f67541bf8f23415146de7"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.11"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "58443b63fb7e465a8a7210828c91c08b92132dff"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.14+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─bd9652fc-f262-11ec-3cf6-71de5b26f0f0
# ╟─8d01fe33-6610-44b1-800a-5fba9224359c
# ╠═38e1904d-b7bb-4868-af51-38688c8f65e6
# ╠═1be1644b-7270-4a8d-a208-8f000df47c5d
# ╠═4cf6bb90-684b-431e-ae72-fdf541550d17
# ╠═c78b8938-2f4e-4b21-afa0-01c6aa03b689
# ╠═64624649-8a70-44fb-bbd5-5c3466ce484a
# ╟─ce747551-460b-4542-b07c-7adc3fdd08cd
# ╠═8a45213f-2429-4035-a782-81e5bf8fba17
# ╟─37940aa2-2876-4d3a-bed4-306aa42d5f42
# ╠═cfc9a88f-4b7e-48e6-800b-830e7d368597
# ╠═370e4f43-3835-45d1-b614-5732c66842f5
# ╟─88564be4-7943-4577-adee-b5e96351912b
# ╠═4a5cfe7c-109f-4bb0-ba17-6a4aba9dc16a
# ╟─287fd231-4ff9-4c2e-a55d-1786a7e41659
# ╠═851b7260-1701-4e92-b2ec-cc9ba37d2e46
# ╟─cdabefa8-77c6-4a56-bef4-c5591330ea59
# ╠═425cb918-4a8b-46e2-810b-17d32c7081f5
# ╟─2908dc16-74ee-4c87-acf8-681b3f2854f2
# ╠═239406a6-0424-4107-a809-2fccb1dfb8b1
# ╠═1fb80604-2385-49d6-b07c-b6b9f467c616
# ╠═5456d0d6-2b69-4a91-9ef2-0ed55361c39f
# ╟─bde7469f-4734-4882-bd3b-2971468dfa3e
# ╠═a42a53d3-8046-4c2a-bbef-4cd4d246b141
# ╟─9812b45d-4973-414e-a71c-22dd81367f88
# ╠═064c6b27-f16d-4fbd-9ec6-67e091ecc9d2
# ╟─4e282e00-0749-4c0a-b44d-ef3cdfbdb10f
# ╠═83ebfa19-b6d6-40db-9cae-bb5e840b4e46
# ╠═36591c12-f038-4ea7-af94-1e60b7aad0ba
# ╠═0af2ae11-cc91-44a2-8fb9-cef16d15fbee
# ╟─df70cfe0-2859-4b22-ba01-ed12db675e27
# ╠═4de0df7b-457d-4a38-a47f-83e536c47473
# ╟─8715a036-c3b3-49bd-af31-d349c0d63bc3
# ╠═b9d102ac-0aab-45f6-9f1d-f23ba590bcc6
# ╟─5e77d309-ca1b-47b2-8c71-a0cb2c96758a
# ╠═96ca485b-db0a-4517-b15b-6cb3232c86e2
# ╟─f5decd43-03ca-40ce-89a0-2cb74fa0fa70
# ╠═5c85e0d5-d0aa-4e32-8ecd-e2688fc51201
# ╠═67bc83af-a427-4380-81fb-dbdb14cd2439
# ╠═58a2aa85-8cce-4986-bfa2-5284afda40e6
# ╠═de790428-9248-4ec3-8bb6-b8fd187483ec
# ╟─95011cf6-e81b-4a48-9fdc-388df2c9ed75
# ╠═c98a4d00-e94d-47aa-aa1b-d91c8422ddee
# ╟─44e657e2-01bf-4fea-a2c3-c9eb900a6214
# ╠═e1f01522-8544-47ca-94ba-0deb9acaba15
# ╟─d3c9d708-f556-4ef6-a34b-c7e461351af4
# ╠═9f2b2856-4256-48dc-8f02-a45dbc898798
# ╟─ad1cec80-6f26-4541-a6af-2363257f609d
# ╠═8529552b-bb06-49bd-937a-0bf21059b144
# ╟─15a5f712-415b-41c4-a3b0-d236075a7aba
# ╠═c641f018-5e61-407c-a651-7c7f7b114449
# ╠═24aa8793-5470-4122-acd7-92a29312689b
# ╠═f322ff0a-49e4-41a2-96d7-e73d5c1f29bd
# ╟─351a621b-c782-45dd-b141-d248eb07e8e8
# ╠═8b6bc20b-6cc6-4450-a8c3-9b58fdbda6d8
# ╟─1e48e79b-7242-473b-a730-0ea7c5313c0f
# ╠═0b80a1b8-b6c6-4663-987d-b4613536056b
# ╠═7fdf6e5c-471e-4170-9a7c-2c027ec02100
# ╠═cd06d809-0a85-485e-9ec7-7698cc748344
# ╠═3a064eac-2e45-4565-a35b-71e3f755c5c0
# ╠═e57c4876-f4c2-40d4-befc-1f6e82ae1496
# ╠═e5cbfd13-7224-453d-b108-a4e54259540d
# ╠═bc0f81b3-3750-426e-8113-712f7490236e
# ╠═f151ec8d-08c4-4e2c-bc2e-555b999d771a
# ╠═71be47ae-1682-4068-a4d7-e36291682efa
# ╟─8d64359d-3ab4-4da0-af92-afc4d25389bf
# ╠═e0675c3d-7e3d-4866-adc0-ade649d3955c
# ╟─42ff9700-8d33-41b0-915a-acb6202bc227
# ╠═82d5d99b-e2b6-46a0-806e-514a019ccb4d
# ╠═d5390ec3-c8c0-4025-92ac-42c762ec361b
# ╠═72a0d2de-aa12-492f-a4f2-22b5eb8e5076
# ╠═bee631cc-5aa4-45e5-9499-38e60bad9418
# ╠═184c8e04-687b-4ea7-90ad-40fb12cd5b13
# ╠═e59580c0-8b15-417a-951c-0e7778a949d9
# ╠═0eabcaa7-4c33-4463-a3bb-c4da1f37b36a
# ╠═9eb41aa2-f534-488a-8f1b-90842353393d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
