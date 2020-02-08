# Installing Slomo

Install Julia (>= 1.3) from the [Julia website](https://julialang.org/).

Download the repository from [GitHub](https://github.com/adwasser/Slomo.jl)

```shell
git clone https://github.com/adwasser/Slomo.jl.git
```

Activate the package and download any dependencies from within Julia

From the shell:

```shell
cd Slomo.jl
julia
```

From the Julia REPL:
	
```julia
julia> using Pkg
julia> Pkg.activate(".")
julia> Pkg.instantiate()
julia> Pkg.resolve()
```

Run the tests:

```julia
julia> using Pkg
julia> Pkg.test("Slomo")
```

(Optional) Build a system image for faster load times.

This requires the [`PackageCompiler.jl` package](https://github.com/JuliaLang/PackageCompiler.jl).

```shell
mkdir build
julia build.jl
```

You can start Julia with the resulting system image by calling

```shell
julia -J /path/to/Slomo.jl/build/sys.so
```

where `/path/to/Slomo.jl` should be replaced by the relevant path to the Slomo installation directory.
