using BEAST
using CompScienceMeshes
using StaticArrays
using Plots
using DelimitedFiles
using LinearAlgebra
using SparseArrays

Ω = CompScienceMeshes.tetmeshsphere(1,0.2)
X = BEAST.nedelecc3d(Ω)
#ttrX = BEAST.ttrace(X,boundary(Ω))
Id = BEAST.Identity()

using Profile
Profile.init()
print(Profile.init(),"\n")
Profile.clear()
@profile assemble(Id,X,X)
Juno.profiletree()
Juno.profiler()
@profiler assemble(Id,X,X)
