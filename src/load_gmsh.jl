using BEAST
using StaticArrays
using DelimitedFiles
using EMwavesBEP
using CompScienceMeshes

d = dirname(pathof(EMwavesBEP))
include(joinpath(d,"gmsh3d.jl"))
Γ = read_gmsh3d_mesh(joinpath(d, "smallfemsphere.msh"))
