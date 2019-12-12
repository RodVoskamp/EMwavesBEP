using BEAST
using CompScienceMeshes
using EMwavesBEP

d = dirname(pathof(EMwavesBEP))
include(joinpath(d,"test_connect.jl"))
print("\n")
include(joinpath(d,"test_ndlcc.jl"))
print("\n")
include(joinpath(d,"tets_ndlcc.jl"))
print("\n")
include(joinpath(d,"test_ndlcdorientation.jl"))
print("\n")
include(joinpath(d,"testcurl.jl"))
print("\n")

print("alles werkt!")
