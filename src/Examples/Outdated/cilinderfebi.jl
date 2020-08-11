using BEAST
using CompScienceMeshes
using StaticArrays
using EMwavesBEP
using Plots
using DelimitedFiles
using LinearAlgebra
using IterativeSolvers

Ω = CompScienceMeshes.tetmeshcilinder(0.1,1,0.1)
#d = dirname(pathof(EMwavesBEP))
#include(joinpath(d,"gmsh3d.jl"))
#Ω = read_gmsh3d_mesh(joinpath(d,"holecilinder.msh"))

Γ = boundary(Ω)
X2 = BEAST.raviartthomas(Γ)
X3 = BEAST.buffachristiansen(Γ)
X = BEAST.nedelecc3d(Ω)
ttrX = BEAST.ttrace(X,Γ)
Y = curl(X)
ttrY = BEAST.ttrace(Y,Γ)

ϵo = 1.0
μo = 1.0
ω = 5.0
κo = ω*sqrt(μo*ϵo)
pfc = im*μo*ω

E = Maxwell3D.planewave(direction=ẑ, polarization=ŷ, wavenumber=κo)
e = (n × E) × n
H = -1/(pfc)*curl(E)
h = (n × H) × n

Id = BEAST.Identity()
N = NCross()
T = Maxwell3D.singlelayer(wavenumber=κo)
K = Maxwell3D.doublelayer(wavenumber=κo)

function κ2(p)
    #ω^2 #PMCHWT
    r = sqrt(p[3]^2+p[2]^2)
    if r > 0.25
        ω^2*(4*(r-0.25)/r)^2
    else
        ω^2*ϵo*μo
    end
end
MP = BEAST.Multiplicative(κ2)

A1 = assemble(Id,Y,Y)
A2 = assemble(MP,X,X)
A3a = assemble(T,ttrX,ttrX)
A3b1 = assemble((K-0.5N),ttrX,X2)
A3b2 = assemble(T,X2,X2)
A3b3 = assemble((K+0.5N),X2,ttrX)

S = (sqrt(ϵo/μo)*A3a+sqrt(ϵo/μo)*(A3b1*inv(A3b2)*A3b3))
Aj = A1-A2-pfc*S

b1 = assemble(e,X3)
b2 = assemble(h,ttrX)

Si1 = assemble(T,ttrX,X2)*inv(Array(assemble(N,X3,X2)))
Si2 = A3b1*inv(A3b2)*assemble((K+0.5N),X2,X2)*inv(Array(assemble(N,X3,X2)))
bj = sqrt(ϵo/μo)*(Si1+Si2)*b1-b2

uj = zeros(length(bj))*im
uj = gmres!(uj, Aj, pfc*bj)
um = -uj/pfc
um = (-S*uj-bj)
ttrY = ttrX
Y = X

using Plots
print("hier")
import PlotlyJS
fcrj, geo = facecurrents(uj,ttrX)
print(extrema(norm.(fcrj)))
display(PlotlyJS.plot(patch(geo, norm.(fcrj))))
fcrm, geo = facecurrents(um,ttrY)
print(extrema(norm.(fcrm)))
display(PlotlyJS.plot(patch(geo, norm.(fcrm))))

t = range(-2,2,length=200)
t2 = range(-2,2,length=200)
pts = [point(0,s,s2) for s in t, s2 in t2]
SLj = potential(BEAST.MWSingleLayerField3D(κo),pts,uj,ttrX)
DLm = potential(BEAST.MWDoubleLayerField3D(κo),pts,um,ttrY)
display(plot(contour(t,t2,norm.(SLj+DLm),camera=(0,90))))
display(plot(heatmap(t,t2,norm.(SLj+DLm),clim=(0, 2),camera=(0,90))))

vals = BEAST.grideval(pts, uj, X)
vals2 = BEAST.grideval(pts, um, Y)
display(plot(surface(t,t2,getindex.(real(vals+vals2),2),camera=(0,90))))
display(plot(surface(t,t2,norm.(vals+vals2),camera=(0,90))))

Φ, Θ = [0.0], range(0,stop=π,length=100)
pts = [point(cos(ϕ)*sin(θ), sin(ϕ)*sin(θ), cos(θ)) for ϕ in Φ for θ in Θ]
#pts = -[point(0.0, sin(θ), cos(θ)) for ϕ in Φ for θ in Θ]
ffd = potential(MWFarField3D(κo), pts, uj, ttrX)
ffd2 = potential(MWFarField3D(κo), pts, um, ttrY)
ff = (ffd2+cross.(pts,ffd))
display(scatter(Θ, real.(norm.(ffd2))))
