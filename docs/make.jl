using Documenter, Microbiome

makedocs(
    sitename = "Microbiome.jl",
    pages = [
        "Home" => "index.md",
        "Microbial Abundances" => "abundances.md",
        "Distances & Dissimilarity" => "distances.md",
        "Contributing" => "contributing.md"
    ],
    authors = "Kevin Bonham, PhD",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true")
)

deploydocs(
    repo = "github.com/BioJulia/Microbiome.jl.git",
<<<<<<< HEAD
    osname = "linux",
    target = "build",
=======
>>>>>>> 20753a3b8f7dfe817bf898b8f4fd13bb4167ee6c
    deps = nothing,
    make = nothing
)
