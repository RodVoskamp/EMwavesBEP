using BEAST
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
tet = simplex(a,b,c,o)
X = BEAST.nedelecc3d(Γ)
X = curl(X)
x = BEAST.refspace(X)

for(q,fc) in enumerate(BEAST.faces(tet))
    Q2 = zeros(scalartype(x),3,4)
    Q = BEAST.ttrace(x,tet,q,fc)
    for i in 1:size(Q,1)
        for j in 1:size(Q,2)
            if Q[i,j] != 0
                Q2[i,j] = 1
            end
        end
    end
    if q == 1
        Q3 = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    elseif q == 2
        Q3 = [0.0 0.0 1.0 0.0; 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0]
    elseif q == 3
        Q3 = [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 0.0 1.0]
    elseif q == 4
        Q3 = [0.0 1.0 0.0 0.0; 1.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0]
    end
    #@assert Q2 == Q3
end
print("untested")
