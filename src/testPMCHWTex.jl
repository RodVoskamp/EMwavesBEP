using CompScienceMeshes
using BEAST

d = dirname(pathof(BEAST))
Γ = readmesh(joinpath(d,"../examples/sphere2.in"))
X = raviartthomas(Γ)

#define test functions

μo, μi, ϵo, ϵi, ω = 1, 1, 1, 1, 1

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

u = s[1:1236,1]
#u = s[1237:2472,1]
using Plots
print("hier")
import PlotlyJS
fcrj, geo = facecurrents(u,X)
print(extrema(norm.(fcrj)))
PlotlyJS.plot(patch(geo, BEAST.norm.(fcrj)))

#u = s[1237:2472,1]
#plotresults = true
#include(joinpath(d,"../examples/utils/postproc.jl"))
#include(joinpath(d,"../examples/utils/plotresults.jl"))
#plot
