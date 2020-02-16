using BEAST
using CompScienceMeshes
using StaticArrays

v = SArray{Tuple{3},Float64,1,3}
f = SArray{Tuple{4},Int,1,4}

b1 = v([4,2,4])
b2 = v([1,13,7])
n1 = v([2,04,12])
n2 = v([0,2,5])
vertice = [b1,b2,n1,n2]

face = zeros(f,1)
face[1] = f([1,2,3,4])

Γ2 = Mesh(vertice,face)
X2 = BEAST.nedelecd3d(Γ2)
t = X2.fns

ref = refspace(X2)
for (i,k) in enumerate(cells(skeleton(Γ2,2)))
    j = setdiff([1,2,3,4],k)
    p = simplex(vertice[k[1]],vertice[k[2]],vertice[k[3]])
    A = volume(p)
    ptch = chart(Γ2,Γ2.faces[1])
    norm2 = carttobary(ptch,0.5*(vertice[k[1]]+vertice[k[2]]))
    ctrd = neighborhood(ptch, norm2)

    #print("  ",1/A,"  ",(dot(t[i][1].coeff*(ref(ctrd)[t[i][1].refid].value),normal(p))),normal(p),"\n")
    @assert abs(1/A-dot(t[i][1].coeff*(ref(ctrd)[t[i][1].refid].value),normal(p)))< 0.001
end
print("test succesvol")
