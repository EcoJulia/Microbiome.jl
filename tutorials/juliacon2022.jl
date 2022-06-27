### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 38e1904d-b7bb-4868-af51-38688c8f65e6
using Downloads

# ╔═╡ cfc9a88f-4b7e-48e6-800b-830e7d368597
using CSV, DataFrames, Microbiome

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
If you have your own data, you should be able to substitute it below
(we'll show you how).
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
ginisimpson!(comm)

# ╔═╡ 24aa8793-5470-4122-acd7-92a29312689b
shannon!(comm)

# ╔═╡ f322ff0a-49e4-41a2-96d7-e73d5c1f29bd
samples(comm)

# ╔═╡ 7fdf6e5c-471e-4170-9a7c-2c027ec02100
md"""
## Working with sample metadata

"""

# ╔═╡ c4df70fd-ea07-4601-846a-b492baed4d60
md"""
## Load data and metadata into Microbiome.jl
"""

# ╔═╡ 5653aecc-4747-4891-bad3-acb800d5774e


# ╔═╡ 42ff9700-8d33-41b0-915a-acb6202bc227
md"""
## Create Principal Coordinate Analysis (PCoA)
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Downloads = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
Microbiome = "3bd8f0ae-a0f2-5238-a5af-e1b399a4940c"

[compat]
CSV = "~0.10.4"
DataFrames = "~1.3.4"
Microbiome = "~0.9.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

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

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "873fb188a4b9d76549b81465b1f75c82aaf59238"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.4"

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

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "9be8be1d8a6f44b96482c8af52238ea7987da3e3"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.45.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

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

[[deps.EcoBase]]
deps = ["RecipesBase"]
git-tree-sha1 = "a4d5b263972e820e780effc2084f92399ba44ee3"
uuid = "a58aae7d-b440-5a11-b283-399458f99aac"
version = "0.1.6"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "129b104185df66e408edd6625d480b7f9e9823a0"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.18"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.Indexing]]
git-tree-sha1 = "ce1566720fd6b19ff3411404d4b977acd4814f9f"
uuid = "313cdc1a-70c2-5d6a-ae34-0150d3930a38"
version = "1.1.1"

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

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

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

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

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

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "0044b23da09b5608b4ecacb4e5e6c6332f833a7e"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

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

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

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

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
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
# ╠═7fdf6e5c-471e-4170-9a7c-2c027ec02100
# ╟─c4df70fd-ea07-4601-846a-b492baed4d60
# ╠═5653aecc-4747-4891-bad3-acb800d5774e
# ╟─42ff9700-8d33-41b0-915a-acb6202bc227
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
