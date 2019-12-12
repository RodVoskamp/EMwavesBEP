using BEAST
using CompScienceMeshes
using StaticArrays

v = SArray{Tuple{3},Float64,1,3}
f = SArray{Tuple{4},Int,1,4}
o = v([0,0,0])
a = v([1,0,0])
b = v([0,2,0])
c = v([0,0,4])
vertice = [a,b,c,o]
face = zeros(f,1)
face[1] = f([1,2,3,4])
Γ = Mesh(vertice,face)
useedg, deledg, delfac = CompScienceMeshes.stripboundedge(Γ)
mesuseedg = useedg.faces
X = BEAST.nedelecc3d(Γ,deledg)
print(curl(X).fns)
X2 = ntrace(curl(X),delfac)
print(X2.fns)
