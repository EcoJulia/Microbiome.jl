using Documenter, Microbiome, Microbiome.Dictionaries

makedocs(
    sitename = "Microbiome.jl",
    pages = [
        "Home" => "index.md",
        "Samples and features" => "samples_features.md",
        "Profiles and Communities" => "profiles.md",
        "Diversity measures" => "diversity.md",
    ],
    authors = "Kevin Bonham, PhD <kbonham@wellesley.edu",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true")
)

deploydocs(
    repo = "github.com/EcoJulia/Microbiome.jl.git",
    push_preview=true,
    devbranch="main"
)
