using BEAST
using CompScienceMeshes
using StaticArrays
using EMwavesBEP
using Plots
using DelimitedFiles
using LinearAlgebra
using IterativeSolvers

Ω = CompScienceMeshes.tetmeshsphere(1,0.2)

#edges = skeleton(Ω, 1)
#faces2 = [sort(c) for c in cells(skeleton(boundary(Ω),2))]
#function is_fint(face)
#    (sort(face) in faces2)
#end

#Γ = submesh(is_fint, boundary(Ω))
Γ = boundary(Ω)
X2 = BEAST.raviartthomas(Γ)
X = BEAST.nedelecc3d(Ω)
ttrX = BEAST.ttrace(X,Γ)
Y = curl(X)

ϵo = 1.0
μo = 1.0
ω = 1.0
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
    #ω^2*(2 - (BEAST.norm(p)^2))
    #ω^2
    #(2 - (BEAST.norm(p)^2))
    1
end
MP = BEAST.Multiplicative(κ2)

A1 = assemble(Id,Y,Y)
A2 = assemble(MP,X,X)
A2b = assemble(Id,X,X)
A3a = assemble(T,ttrX,ttrX)
A3b1 = assemble((K+0.5N),ttrX,X2)
#A3b1 = assemble((K-0.5N),ttrX,X2)
A3b2 = assemble(T,X2,X2)
A3b3 = assemble((K-0.5N),X2,ttrX)
#A3b3 = assemble((K+0.5N),X2,ttrX)

#S = (ϵo/μo*A3a+ϵo/μo*(A3b1*(A3b2)*A3b3))
S = (ϵo/μo*A3a+ϵo/μo*(A3b1*inv(A3b2)*A3b3))

Aj = A1-A2-pfc*S*A2b
#Aj = A1-A2+pfc*S*A2b
#Am = A1-A2+(-pfc*-S*A2b)

b1 = assemble(e,ttrX)
b2 = assemble(h,ttrX)
bj = -pfc*S*b1+pfc*b2
#bj = pfc*S*b1+pfc*b2
#bm = (-pfc)*(-S)*b1-pfc*b2

uj = zeros(length(bj))*im
#um = zeros(length(bm))*im
#uj = Aj\bj
#um = Am\bm
uj = gmres!(uj, Aj, bj)
#um = gmres!(um, Am, bm)
um = S*(uj-b1)+b2

using Plots
print("hier")
import PlotlyJS
fcrj, geo = facecurrents(uj,ttrX)
print(extrema(norm.(fcrj)))
display(PlotlyJS.plot(patch(geo, norm.(fcrj))))
fcrj, geo = facecurrents(um,ttrX)
print(extrema(norm.(fcrj)))
display(PlotlyJS.plot(patch(geo, norm.(fcrj))))

t = range(-2,2,length=200)
t2 = range(-2,2,length=200)
pts = [point(s,0,s2) for s in t, s2 in t2]
SLj = potential(BEAST.MWSingleLayerField3D(κo),pts,uj,ttrX)
DLm = potential(BEAST.MWDoubleLayerField3D(κo),pts,um,ttrX)
display(plot(contour(t,t2,norm.(SLj+DLm),camera=(0,90))))
display(plot(heatmap(t,t2,norm.(SLj+DLm),clim=(-4*ω, 4*ω),camera=(0,90))))
display(plot(heatmap(t,t2,getindex.(real(SLj+DLm),1),clim=(-4, 4),camera=(0,90))))
display(plot(heatmap(t,t2,getindex.(real(SLj),1),clim=(-2*ω, 2*ω),camera=(0,90))))

const CSM = CompScienceMeshes
hemi = submesh(tet -> cartesian(CSM.center(chart(Ω,tet)))[2] < 0, Ω)
bnd_hemi = boundary(hemi)

Xhemi = BEAST.restrict(X, hemi)
tXhemi = BEAST.ttrace(Xhemi, bnd_hemi)

fcr, geo = facecurrents(uj, tXhemi)
PlotlyJS.plot(patch(geo, norm.(fcr)))

#small test
S = (ϵo/μo*A3a+ϵo/μo*(A3b1*inv(A3b2)*A3b3))
S2 = norm.(S[1,:])
S3 = zeros(1134,1134)
l = 0
j = []
for i = 1:3117
    if S2[i] ≈ 0
        global l += 1
        append!(j,i)
    end
end
S = S[setdiff(1:end, j), :]
S = S[:, setdiff(1:end, j)]
print(3117-1134," ",l)
