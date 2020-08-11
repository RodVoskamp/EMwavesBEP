using CompScienceMeshes
using BEAST
using LinearAlgebra
using EMwavesBEP
using StaticArrays
using DelimitedFiles

#d = dirname(pathof(BEAST))
#Γ = readmesh(joinpath(d,"../examples/sphere2.in"))
Ω = CompScienceMeshes.tetmeshsphere(1,0.12)
#Ω = CompScienceMeshes.tetmeshcilinder(1,1,0.3)
#d = dirname(pathof(EMwavesBEP))
#include(joinpath(d,"gmsh3d.jl"))
#Ω = read_gmsh3d_mesh(joinpath(d,"holecilinder.msh"))

Γ = boundary(Ω)
X = raviartthomas(Γ)

#define test functions

μo, μi, ϵo, ϵi, ω = 1, 1, 1, 1, 5

κ, κin = ω*sqrt(ϵo*μo), ω*sqrt(ϵi*μi)
η = sqrt(μi/ϵi)

#define variables

t = Maxwell3D.singlelayer(wavenumber=κ)
tin = Maxwell3D.singlelayer(wavenumber=κin)
v = Maxwell3D.doublelayer(wavenumber=κ)
vin = Maxwell3D.doublelayer(wavenumber=κin)

#matrix operators

E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
H = -1/(im*ω*μo)*curl(E)
e = (n × E) × n
h = (n × H) × n

# f and g

@hilbertspace k l
@hilbertspace j i

#it will be an [A,B;C,D][I;K]=[f;g], so k+l by j+i

A = (t+η*tin)
C = (v+vin)
B = (-1*C)
D = (t+1/η*tin)

pmchwt = @discretise A[k,j]+B[k,i]+C[l,j]+D[l,i]==e[k]+h[l] k∈X l∈X j∈X i∈X
s = solve(pmchwt)

#s = [I;K]
u2 = s[1:Int(length(s)/2),1]
u = s[Int(length(s)/2)+1:length(s),1]
using Plots
print("hier")
import PlotlyJS
fcrj, geo = facecurrents(u,X)
a,b,c = assemblydata(X)
pts = CompScienceMeshes.center.(a)
EE = real.(E.(cartesian.(pts)))
print(extrema(norm.(fcrj)))
display(PlotlyJS.plot(patch(geo, BEAST.norm.(fcrj))))
fcrj, geo = facecurrents(u2,X)
print(extrema(norm.(fcrj)))
display(PlotlyJS.plot(patch(geo, BEAST.norm.(fcrj))))

Φ, Θ = [0.0], range(0,stop=π,length=100)
pts = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for ϕ in Φ for θ in Θ]
ffd = potential(MWFarField3D(κ), pts, u2, X)
display(scatter(Θ, real.(norm.(ffd))))
