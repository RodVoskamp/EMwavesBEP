using BEAST
using CompScienceMeshes
using StaticArrays
using EMwavesBEP

d = dirname(pathof(BEAST))
include(joinpath(d,"../examples/efie.jl"))

d2 = dirname(pathof(EMwavesBEP))
include(joinpath(d2,"load_gmsh.jl"))

X = BEAST.nedelecd3d(Γ)
