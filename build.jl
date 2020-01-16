using PackageCompiler

target = joinpath("build", "sys.so")

if ispath(target)
    mv(target, target * ".bak"; force = true)
end

imgfile = compile_package("DifferentialEquations", "Cosmology",
                          "HypergeometricFunctions", "Interpolations", "Roots")
filename = splitdir(imgfile)[2]
@assert filename == splitdir(target)[2]
mv(imgfile, target)
