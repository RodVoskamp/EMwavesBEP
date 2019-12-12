using BEAST
using CompScienceMeshes

Γ = CompScienceMeshes.tetmeshsphere(1,0.4);
totedg = size(skeleton(Γ,1).faces,1)
useedg, deledg, delfac = CompScienceMeshes.stripboundedge(Γ)
totuseedg, totdeledg, totdelfac = size(useedg.faces,1), size(deledg.faces,1), size(delfac.faces,1)
@assert totuseedg+totdeledg==totedg
@assert size(setdiff(useedg.vertices,deledg.vertices),1) == 0
@assert size(setdiff(useedg.vertices,delfac.vertices),1) == 0

C1 = CompScienceMeshes.connectivity(useedg,delfac)
C2 = CompScienceMeshes.connectivity(deledg,delfac)
@assert size(CompScienceMeshes.rowvals(C1),1) == 0
@assert size(CompScienceMeshes.rowvals(C2),1) == 3*size(delfac.faces,1)
@assert size(CompScienceMeshes.rowvals(C2),1)/2 == size(skeleton(delfac,1).faces,1)

print("test succesvol")
