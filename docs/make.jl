using Documenter, Microbiome

makedocs()

deploydocs(
    repo = "github.com/BioJulia/Microbiome.jl.git",
    julia = "0.6",
    osname = "linux",
    target = "build",
    deps = nothing,
    make = nothing
)
