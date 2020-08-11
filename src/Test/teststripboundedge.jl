using BEAST
using CompScienceMeshes
using EMwavesBEP
using StaticArrays
using PyPlot
using Plots
using DelimitedFiles

pyplot()

d = dirname(pathof(EMwavesBEP))
#include(joinpath(d,"stripboundedge.jl"))
include(joinpath(d,"gmsh3d.jl"))

pygui(true)
fig = figure()
ax = fig.gca(projection="3d")

Γ = read_gmsh3d_mesh(joinpath(d,"smallcube.msh"))
Γ = CompScienceMeshes.tetmeshsphere(3,0.5)
vert = Γ.vertices
usededges, dummy, dummy2 = CompScienceMeshes.stripboundedge(Γ)
usededges = usededges.faces

l = length(usededges)
i = 0
while i < l
    edge = usededges[i+1]
    if edge[1] != edge[2]
        a = vert[edge[1]]
        b = vert[edge[2]] - a
        ax.quiver(a[1],a[2],a[3], b[1],b[2],b[3],color="red")
    end
    global i += 1
end

fig.canvas.draw()
