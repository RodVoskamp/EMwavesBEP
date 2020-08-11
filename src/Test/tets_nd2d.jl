using BEAST
using CompScienceMeshes
using StaticArrays

v = SArray{Tuple{3},Float64,1,3}
f = SArray{Tuple{3},Int,1,3}

b1 = v([-1,0,0])
b2 = v([1,0,0])
n1 = v([0,1,0])
n2 = v([0,-1,0])
vertice = [b1,b2,n1,n2]

face = zeros(f,2)
for i = 1:length(vertice)-3
    face[i] = f([1,2,i+2])
end
face[2] = f([4,1,2])

Γ2 = Mesh(vertice,face)
X2 = BEAST.nedelec(Γ2)

a = length.(X2.fns)
k = 0.0
shapepos = 0
while k < 1.5
    global shapepos = shapepos+1
    global k = a[shapepos]
end
t = X2.fns[shapepos]

ref = refspace(X2)
l = BEAST.norm(b1-b2)
for (i,k) in enumerate(cells(Γ2))
    #print(k)
    ptch = chart(Γ2,k)
    tang2 = carttobary(ptch,0.5*(b2+b1))
    ctrd = neighborhood(ptch, tang2)

    #print(abs(1/l-dot(t[i].coeff*(ref(ctrd)[t[i].refid][1]),BEAST.normalize(b2-b1)))< 0.001)
    @assert abs(1/l-dot(t[i].coeff*(ref(ctrd)[t[i].refid][1]),BEAST.normalize(b2-b1)))< 0.001
end
print("test succesvol")
