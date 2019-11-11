using BEAST
using CompScienceMeshes
using StaticArrays

'''
#Makes the mesh: b1->b2 is the edge every tetrahedron is centered on
'''

v = SArray{Tuple{3},Float64,1,3}
f = SArray{Tuple{4},Int,1,4}

b1 = v([0,0,0])
b2 = v([0,0,1])
n1 = v([0,1,0.5])
n2 = v([1,0,0.5])
n3 = v([0,-1,0.5])
n4 = v([-1,0,0.5])
vertice = [b1,b2,n1,n2,n3,n4]

face = zeros(f,4)
for i = 1:length(vertice)-3
    face[i] = f([1,2,i+2,i+3])
end
face[4] = f([1,3,length(vertice),2])
#THis edge should have shape{4,2,1}, ans this is what happens.

Γ2 = Mesh(vertice,face)
X2 = BEAST.nedelecc3d(Γ2)

'''
#The next part gives the shape of our edge b1->b2 (t)
'''

a = length.(X2.fns)
k = 0.0
shapepos = 0
while k < 3.5
    global shapepos = shapepos+1
    global k = a[shapepos]
end
t = X2.fns[shapepos]

'''
#Try to test if our previous code is correct
'''
function refid(a)
    if a < 4
        z = point(0,0,0)
    elseif a ≈ 6
        z = point(0,1,0)
    else
        z = point(1,0,0)
    end
end

ref = refspace(X2)
l = BEAST.norm(b1-b2)
for (i,k) in enumerate(cells(Γ2))
    print(k)
    ptch = chart(Γ2,k)
    ctrd = neighborhood(ptch, 0.4*(b2-b1))
    tu = tangents(ctrd,1)
    tv = tangents(ctrd,2)
    tw = tangents(ctrd,3)
    B = [tv-tu tw-tu -tu]
    tang = inv(B)*(0.4*(b2-b1))+refid(t[i].refid)
    tang2 = carttobary(ptch,0.4*(b2-b1))
    print(tang, tang2)
    ctrd = neighborhood(ptch, tang)

    print(dot(t[i].coeff*(ref(ctrd)[t[i].refid].value),BEAST.normalize(b2-b1)))
    print("=", 1/BEAST.norm(b2-b1))
    print("\n")
end
