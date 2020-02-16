using CompScienceMeshes
using BEAST

tetrs = CompScienceMeshes.tetmeshsphere(1.0, 0.40)
#tetrs = CompScienceMeshes.tetmeshcuboid(1.0, 2.0, 3.0, 1)
@show numcells(tetrs)

bndry = boundary(tetrs)
edges = skeleton(tetrs, 1)
faces = skeleton(tetrs, 2)

bndry_edges = [sort(c) for c in cells(skeleton(bndry, 1))]
bndry_faces = [sort(c) for c in cells(skeleton(bndry, 2))]
function is_interior(edge)
    !(sort(edge) in bndry_edges)
end
function is_fint(face)
    (sort(face) in bndry_faces)
end

interior_edges = submesh(is_interior, edges)
delfac = submesh(is_fint, faces)
@assert numcells(interior_edges) + numcells(skeleton(bndry,1)) == numcells(edges)

X = BEAST.nedelecc3d(tetrs, interior_edges)
@assert numfunctions(X) == numcells(interior_edges)

Id = BEAST.Identity()
Y = curl(X)
A1 = assemble(Id, Y, Y)
A2 = assemble(Id, X, X)
A = A1 - 10*A2

using LinearAlgebra
f = BEAST.ScalarTrace(x -> point(1,0,0) * exp(-norm(x)^2/4))
b = assemble(f, X)
print("hier")
u = A \ b

using Plots

#pts = [point(0,t,0) for t in range(-2,2,length=200)];
#vals = BEAST.grideval(pts, u, X)
#plot(getindex.(real(vals),1), m=2)
print("hier")
import PlotlyJS
ttrY = BEAST.ttrace(Y,delfac)
fcrj, geo = facecurrents(u,ttrY)
PlotlyJS.plot(patch(geo, norm.(fcrj)))
