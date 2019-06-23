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
    osname = "linux",
    target = "build",
    deps = nothing,
    make = nothing
)
