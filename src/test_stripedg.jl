using BEAST
using EMwavesBEP
using CompScienceMeshes
using StaticArrays

v = SArray{Tuple{3},Float64,1,3}
f = SArray{Tuple{4},Int,1,4}

o = v([0,0,0])
a = v([1,0,0])
b = v([0,2,0])
c = v([0,0,4])
face = zeros(f,1)
face[1] = f([1,2,3,4])
Γ = Mesh([a,b,c,o], face)

useedg, deledg, delfac = CompScienceMeshes.stripboundedge(Γ)
@assert size(useedg.faces,1) == 0
@assert size(deledg.faces,1) == 6
@assert size(delfac.faces,1) == 4


b1 = v([-1,0,0])
b2 = v([1,0,0])
n1 = v([0,1,0])
n2 = v([0,0,-1])
n3 = v([0,-1,0])
n4 = v([0,0,1])
vertice = [b1,b2,n1,n2,n3,n4]

face = zeros(f,4)
for i = 1:length(vertice)-3
    face[i] = f([1,2,i+2,i+3])
end
face[4] = f([3,1,2,length(vertice)])
Γ2 = Mesh(vertice,face)

useedg, deledg, delfac = CompScienceMeshes.stripboundedge(Γ2)
@assert size(useedg.faces,1) == 1
@assert size(deledg.faces,1) == 12
@assert size(delfac.faces,1) == 8

print("test succesvol")
