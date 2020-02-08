push!(LOAD_PATH,"../src/")

using Documenter, Slomo

makedocs(
    sitename = "Slomo.jl",
    modules = [Slomo],
    doctest = false,
    pages = [
	"Getting started" => "index.md",
	"Installation" => "install.md",
        "Modules" => [
	    "modules/models.md",
	    "modules/tracers.md",
	    "modules/halos.md",
	    "modules/anisotropy.md",
	    "modules/jeans.md"
        ]
    ]
)
