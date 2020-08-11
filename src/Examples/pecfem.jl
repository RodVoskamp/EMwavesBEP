using CompScienceMeshes
using BEAST

tetrs = CompScienceMeshes.tetmeshsphere(1.0, 0.10)
#tetrs = CompScienceMeshes.tetmeshcuboid(1.0, 2.0, 3.0, 0.3)
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
A = A1 - 9*A2

using LinearAlgebra
f = BEAST.ScalarTrace(x -> point(1,0,0) * exp(-norm(x)^2/4))
b = assemble(f, X)
print("hier")
u = A \ b

using Plots

pts = [point(0,0,t) for t in range(-1.5,1.5,length=20000)];
zx = range(-1.5,1.5,length=20000)
vals = BEAST.grideval(pts, u, X)
plot(zx,getindex.(real(vals),1),label = ["x"])
plot!(zx,getindex.(real(vals),2),label = ["x","y"])
display(plot!(zx,getindex.(real(vals),3),label = ["x","y","z"], xlabel="z", title = "value for the x, y and z direction in (0,0,z)'"))

print("hier")
import PlotlyJS
ttrY = BEAST.ttrace(Y,delfac)
fcrj, geo = facecurrents(u,ttrY)
mini,maxi = extrema(norm.(fcrj))
display(PlotlyJS.plot(patch(geo, norm.(fcrj),(mini,maxi))))

const CSM = CompScienceMeshes
hemi = submesh(tet -> cartesian(CSM.center(chart(tetrs,tet)))[2] < 0, tetrs)
bnd_hemi = boundary(hemi)

Xhemi = BEAST.restrict(X, hemi)
tXhemi = BEAST.ttrace(Xhemi, bnd_hemi)
fcr, geo = facecurrents(u, tXhemi)
display(PlotlyJS.plot(patch(geo, norm.(fcr))))

Yhemi = BEAST.restrict(Y, hemi)
tYhemi = BEAST.ttrace(Yhemi, bnd_hemi)
fcr, geo = facecurrents(u, tYhemi)
display(PlotlyJS.plot(patch(geo, norm.(fcr),(mini,maxi))))

timesteps = 100
t = range(-2,2,length=200)
t2 = range(-2,2,length=200)
pts = [point(s,0,s2) for s in t, s2 in t2]
vals = BEAST.grideval(pts, u, X)
tijd = range(0,stop=2π,length=timesteps)
vt = getindex.(vals,1)

function valt(tijd)
    real(exp.(im*tijd)*vt)
end

vv = zeros(0)
for time in tijd
    global vv
    append!(vv,valt(time))
end
0 = 9
anim = @animate for i ∈ 1:timesteps
    plot(surface(t,t2,vv[length(vt)*(i-1)+1:length(vt)*i],camera=(0,90),clim=(-2, 2)))
end
gif(anim, "anim_k9fem.gif", fps = 15)
