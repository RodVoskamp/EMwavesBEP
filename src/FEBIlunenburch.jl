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
faces2 = [sort(c) for c in cells(skeleton(boundary(Ω),2))]
function is_fint(face)
    (sort(face) in faces2)
end
Γ = submesh(is_fint, boundary(Ω))

X2 = BEAST.raviartthomas(Γ)
X = BEAST.nedelecc3d(Ω)
ttrX = BEAST.ttrace(X,Γ)
Y = curl(X)

ϵo = 1.0
μo = 1.0
ω = 3.0
κo = ω*sqrt(μo*ϵo)
pfc = im*μo*ω

E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κo)
e = (n × E) × n
H = -1/(pfc)*curl(E)
h = (n × H) × n

Id = BEAST.Identity()
N = NCross()
T = Maxwell3D.singlelayer(wavenumber=κo)
K = Maxwell3D.doublelayer(wavenumber=κo)

function κ2(p)
    ω^2*(2 - (BEAST.norm(p)^2))
    #1
end
MP = BEAST.Multiplicative(κ2)

A1 = assemble(Id,Y,Y)
A2 = assemble(MP,X,X)
A2b = assemble(Id,X,X)
A3a = assemble(T,ttrX,ttrX)
A3b1 = assemble((K+0.5N),ttrX,X2)
#A3b1 = assemble((K-0.5N),BEAST.ttrace(X,Γ),X2)
A3b2 = assemble(T,X2,X2)
A3b3 = assemble((K-0.5N),X2,ttrX)
#A3b1 = assemble((K+0.5N),BEAST.ttrace(X,Γ),X2)
#S = (ϵo/μo*A3a+ϵo/μo*(A3b1*A3b2*A3b3))
S = (ϵo/μo*A3a+ϵo/μo*(A3b1*inv(A3b2)*A3b3))

#AA = A1-A2-pfc*S*A2
AA = A1-A2+pfc*S*A2b

b1 = assemble(e,ttrX)
b2 = assemble(h,ttrX)
#b = -pfc*S*b1+pfc*b2
b = pfc*S*b1+pfc*b2

u = AA\b

using Plots
print("hier")
import PlotlyJS
fcrj, geo = facecurrents(u,ttrX)
print(extrema(norm.(fcrj)))
PlotlyJS.plot(patch(geo, norm.(fcrj)))

#scatter(norm.(S[110,:]),norm.(S[:,110]))

#pts = [point(0,0,t) for t in range(-2,2,length=1200)];
#vals = BEAST.grideval(pts, u, X)
#plot(getindex.(real(vals),1), m=2)
