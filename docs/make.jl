push!(LOAD_PATH,"../src/")

using Documenter, Slomo

makedocs(
    sitename = "Slomo.jl",
    modules = [Slomo],
    doctest = false,
)
