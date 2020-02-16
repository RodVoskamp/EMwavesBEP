using BEAST
using CompScienceMeshes
using StaticArrays
using EMwavesBEP
using Plots
using DelimitedFiles
using LinearAlgebra

Ω = CompScienceMeshes.tetmeshsphere(1,0.2)
#useedg, deledg, Γ = CompScienceMeshes.stripboundedge(Ω)

Γ = boundary(Ω)
print(CompScienceMeshes.isoriented(Γ))

edges = skeleton(Ω, 1)
faces = skeleton(Ω, 2)
edges2 = [sort(c) for c in cells(edges)]
faces2 = [sort(c) for c in cells(skeleton(boundary(Ω),2))]
function is_fint(face)
    (sort(face) in faces2)
end
Γ = submesh(is_fint, boundary(Ω))
X2 = BEAST.raviartthomas(Γ)
ttrX = BEAST.ttrace(X,Γ)
bc = BEAST.buffachristiansen(Γ)
X = BEAST.nedelecc3d(Ω)
Y = curl(X)

ϵi, ϵo = 1.0, 1.0
μi, μo = 01.0, 1.0
ω = 1.0
κi, κo = ω*sqrt(μi*ϵi), ω*sqrt(μo*ϵo)
pfc = im*μo*ω
pfci = im*μi*ω

E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κo)
e = (n × E) × n
H = -1/(pfc)*curl(E)
h = (n × H) × n

Id = BEAST.Identity()
N = NCross()
T = Maxwell3D.singlelayer(wavenumber=κo)
K = Maxwell3D.doublelayer(wavenumber=κo)

A1 = assemble(Id,Y,Y)
A2 = assemble(Id,X,X)
A3a = assemble(T,ttrX,ttrX)
A3b1 = assemble((K+0.5N),ttrX,bc)
A3b2 = assemble(T,bc,bc)
A3b3 = assemble((K-0.5N),bc,ttrX)
S = (ϵo/μo*A3a+ϵo/μo*(A3b1*inv(A3b2)*A3b3))

AA = A1-κi^2*A2+pfc*S*A2

b1 = assemble(e,ttrX)
b2 = assemble(h,ttrX)
b = pfc*S*b1+pfc*b2

u = AA\b

using Plots
print("hier")
import PlotlyJS
fcrj, geo = facecurrents(u,ttrX)
PlotlyJS.plot(patch(geo, norm.(fcrj)))

#pts = [point(0,t,0) for t in range(-2,2,length=200)];
#vals = BEAST.grideval(pts, u, X)
#plot(getindex.(real(vals),1), m=2)
