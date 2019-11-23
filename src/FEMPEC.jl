using BEAST
using StaticArrays
using DelimitedFiles
using EMwavesBEP
using CompScienceMeshes

d = dirname(pathof(EMwavesBEP))

include(joinpath(d,"gmsh3d.jl"))
include(joinpath(d,"stripboundedge.jl"))
Γ = read_gmsh3d_mesh(joinpath(d,"smallfemsphere.msh"))
useedg, deledg = stripboundedge(Γ)
mesuseedg = Mesh(Γ.vertices,useedg)
X = BEAST.nedelecc3d(Γ,mesuseedg)

I = BEAST.Identity()
assemble(I,curl(X),curl(X))
#assemble(I,X,X)

function f(r) return 0 end

@hilbertspace e
@hilbertspace e2
@discretise I[curl(e2),curl(e)]-I[e2,e] == f[e2] e∈X e2∈X
