using BEAST
using StaticArrays
using DelimitedFiles
using EMwavesBEP
using CompScienceMeshes

Γ = CompScienceMeshes.tetmeshsphere(1,0.4)
#d = dirname(pathof(EMwavesBEP))
#include(joinpath(d,"load_gmsh.jl"))
useedg, deledg, delfac = CompScienceMeshes.stripboundedge(Γ)
X = BEAST.nedelecc3d(Γ,useedg)

using LinearAlgebra
Id = BEAST.Identity()
f(x) = point(0,0,1)*exp(-norm(x)^2)
F =  BEAST.SourceField(f)

@hilbertspace e
@hilbertspace e2
Eq = @varform Id[curl(e2),curl(e)]-1*Id[e2,e] == F[e2]
eq = @discretise Eq e∈X e2∈X
u = solve(eq)
@assert size(useedg.faces,1) == size(u,1)

import PlotlyJS
Y = curl(X)
@assert size(Y.pos,1) == size(u,1)

ttrY = BEAST.ttrace(Y,delfac)
fcrj, _ = facecurrents(u,ttrY)
PlotlyJS.plot(patch(ttrY.geo, norm.(fcrj)))
PlotlyJS.plot(patch(delfac, norm.(fcrj)))

#ntrY = BEAST.ntrace(Y,delfac)
#@assert size(ntrY.pos,1) == size(u,1)
#fcrj, _ = facecurrents(u,ntrY)
#@assert size(fcrj,1) == size(skeleton(Γ,2).faces,1)
#PlotlyJS.plot(patch(skeleton(Γ,2), norm.(fcrj)))

#using Plots
#pts = [point(t,0,0) for t in range(-2,2,length=200)];
#vals = BEAST.grideval(pts, u, X)
#plot(getindex.(norm.(vals),1))
