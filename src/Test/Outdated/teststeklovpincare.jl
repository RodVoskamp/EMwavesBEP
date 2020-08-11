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
Γ = submesh(is_fint, boundary(Ω))
X2 = BEAST.raviartthomas(Γ)
X = BEAST.nedelecc3d(Ω)
Y = curl(X)

ϵi, ϵo = 1.0, 1.0
μi, μo = 1.0, 1.0
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
T = Maxwell3D.singlelayer(wavenumber=κi)
K = Maxwell3D.doublelayer(wavenumber=κi)

#A3a = assemble(T,BEAST.ttrace(X,Γ),BEAST.ttrace(X,Γ))
#A3b1 = assemble((K+0.5N),BEAST.ttrace(X,Γ),X2)
#A3b2 = assemble(T,X2,X2)
#A3b3 = assemble((K-0.5N),X2,BEAST.ttrace(X,Γ))
A3a = assemble(T,X2,X2)
A3b1 = assemble((K+0.5N),X2,X2)
A3b2 = assemble(T,X2,X2)
A3b3 = assemble((K-0.5N),X2,X2)

#ji = assemble(h,BEAST.ttrace(X,Γ))
#mi = assemble(e,BEAST.ttrace(X,Γ))
ji = assemble(h,X2)
mi = assemble(e,X2)

So = (ϵo/μo*A3a+ϵo/μo*(A3b1*inv(A3b2)*A3b3))
Si = (ϵi/μi*A3a+ϵi/μi*(A3b1*inv(A3b2)*A3b3))

η = sqrt(μi/ϵi)
t = Maxwell3D.singlelayer(wavenumber=κo)
tin = Maxwell3D.singlelayer(wavenumber=κi)
v = Maxwell3D.doublelayer(wavenumber=κo)
vin = Maxwell3D.doublelayer(wavenumber=κi)
@hilbertspace k l
@hilbertspace j i
A = (t+η*tin)
C = (v+vin)
B = (-1*C)
D = (t+1/η*tin)
pmchwt = @discretise A[k,j]+B[k,i]+C[l,j]+D[l,i]==e[k]+h[l] k∈X2 l∈X2 j∈X2 i∈X2
s = solve(pmchwt)
u = s[1:Int(length(s)/2),1]
u = s[Int(length(s)/2)+1:length(s),1]

using Plots
print("hier")
import PlotlyJS
fcrj, geo = facecurrents(u,X2)
PlotlyJS.plot(patch(geo, BEAST.norm.(fcrj)))

#@assert So*ji-mi == (So+Si)*u
print("\n",((So+Si)*u)[100],"\n",(So*ji-mi)[100],"\n",-1*mi[100],"\n",(So*ji)[100])
