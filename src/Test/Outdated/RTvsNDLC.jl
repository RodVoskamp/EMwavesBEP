using BEAST
using CompScienceMeshes

Γ = CompScienceMeshes.tetmeshsphere(1,0.3)
useedg, deledg, delfac = CompScienceMeshes.stripboundedge(Γ)

ΓRT = delfac
#ΓRT = CompScienceMeshes.meshsphere(1,0.3)
XRT = raviartthomas(ΓRT)
κ = 1.0
t = Maxwell3D.singlelayer(wavenumber=κ)
E = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=κ)
e = (n × E) × n
@hilbertspace j
@hilbertspace k
efie = @discretise t[k,j]==e[k]  j∈XRT k∈XRT
urt = solve(efie)
import PlotlyJS
using LinearAlgebra
fcrj, _ = facecurrents(urt,XRT)
PlotlyJS.plot(patch(ΓRT, norm.(fcrj)))
