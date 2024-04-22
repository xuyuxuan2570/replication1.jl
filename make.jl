push!(LOAD_PATH, "../src/") 

using Documenter
include("../src/parameter.jl") 

makedocs(sitename="My Documentation")  

