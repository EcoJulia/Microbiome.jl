using Documenter, Microbiome

makedocs(
    format = :html,
    sitename = "Microbiome.jl",
    pages = [
        "Home" => "index.md",
        "Microbial Abundances" => "abundances.md",
        "Distances & Dissimilarity" => "distances.md",
        "Contributing" => "contributing.md"
    ],
    authors = "Kevin Bonham, PhD"
)

deploydocs(
    repo = "github.com/BioJulia/Microbiome.jl.git",
    julia = "0.6",
    osname = "linux",
    target = "build",
    deps = nothing,
    make = nothing
)
