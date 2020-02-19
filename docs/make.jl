src_path = joinpath(splitdir(dirname(@__FILE__))[1], "src")
push!(LOAD_PATH, src_path)

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

deploydocs(
    repo = "github.com/adwasser/Slomo.jl.git"
)

