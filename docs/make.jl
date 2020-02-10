push!(LOAD_PATH,"../src/")

using Documenter, Slomo

makedocs(
    sitename = "Slomo.jl",
    modules = [Slomo],
    doctest = false,
    pages = [
	"Home" => "index.md",
	"Installation" => "install.md",
        "Getting started" => "start.md",
        "Examples" => "examples.md",
        "Public API" => "public_api.md",
        "Private API" => "private_api.md"
    ]
)
