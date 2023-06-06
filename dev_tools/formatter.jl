using Pkg
Pkg.activate(".")
using JuliaFormatter

format("../src")
format("../test")

println("Successfully formatted all files in thie repo.")
