using Documenter, Microbiome, Microbiome.Dictionaries

makedocs(
    sitename = "Microbiome.jl",
    warnonly = [:missing_docs],
    pages = [
        "Home" => "index.md",
        "Samples and features" => "samples_features.md",
        "Profiles and Communities" => "profiles.md",
        "Diversity measures" => "diversity.md",
    ],
    authors = "Kevin Bonham, PhD <kbonham@wellesley.edu>",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        edit_link = "main",
        canonical = "https://docs.ecojulia.org/Microbiome.jl/stable/",
    ),
)

deploydocs(
    repo = "github.com/EcoJulia/Microbiome.jl.git",
    push_preview=true,
    devbranch="main"
)
