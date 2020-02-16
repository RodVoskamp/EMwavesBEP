using BEAST
using CompScienceMeshes
using StaticArrays
using EMwavesBEP
using Plots
using DelimitedFiles
using LinearAlgebra

Ω = CompScienceMeshes.tetmeshsphere(1,0.2)
#useedg, deledg, Γ = CompScienceMeshes.stripboundedge(Ω)

edges = skeleton(Ω, 1)
faces = skeleton(Ω, 2)
edges2 = [sort(c) for c in cells(edges)]
faces2 = [sort(c) for c in cells(skeleton(boundary(Ω),2))]
function is_interior(edge)
    (sort(edge) in edges2)
end
function is_fint(face)
    (sort(face) in faces2)
end
useedg = submesh(is_interior, edges)
Γ = submesh(is_fint, faces)
X2 = BEAST.raviartthomas(Γ)
X = BEAST.nedelecc3d(Ω, useedg)
Y = curl(X)

ϵi, ϵo = 1.0, 1.0
μi, μo = 1.0, 1.0
ω = 1.0
κi, κo = ω*sqrt(μi*ϵi), ω*sqrt(μo*ϵo)
pfc = im*μo*ω

h = BEAST.ScalarTrace(x -> point(1,0,0) * exp(-norm(x)^2/4))

Id = BEAST.Identity()
N = NCross()
T = Maxwell3D.singlelayer(wavenumber=κi)
K = Maxwell3D.doublelayer(wavenumber=κi)

A1 = assemble(Id,Y,Y)
A2 = assemble(Id,X,X)
A3a = assemble(T,BEAST.ttrace(X,Γ),BEAST.ttrace(X,Γ))
A3b1 = assemble((K+0.5N),BEAST.ttrace(X,Γ),X2)
A3b2 = assemble(T,X2,X2)
A3b3 = assemble((K-0.5N),X2,BEAST.ttrace(X,Γ))
S = (ϵo/μo*A3a+ϵo/μo*(A3b1*A3b2*A3b3))
Ax = assemble(Id,BEAST.ttrace(X,Γ),BEAST.ttrace(X,Γ))

A = A1-κi^2*A2-pfc*S
A = A1-κi^2*A2+pfc*S*Ax
A = A1-κi^2*A2+pfc*S*A2

b2 = assemble(h,X)
b = b2

u = A\b

using Plots
print("hier")
import PlotlyJS
ttrX = BEAST.ttrace(X,Γ)
fcrj, geo = facecurrents(u,ttrX)
PlotlyJS.plot(patch(geo, norm.(fcrj)))

#pts = [point(0,0,t) for t in range(-2,2,length=200)];
#vals = BEAST.grideval(pts, u, X)
#plot(getindex.(real(vals),1), m=2)
