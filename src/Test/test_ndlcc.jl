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

edgee = skeleton(Γ,1).faces

function tempor(edg)
    a = min(edg[1],edg[2])
    b = max(edg[1],edg[2])
    if a == 1 && b == 4
        return 3
    elseif a == 2 && b == 4
        return 5
    elseif a == 3 && b == 4
        return 6
    elseif a == 1 && b == 2
        return 1
    elseif a == 1 && b == 3
        return 2
    elseif a == 2 && b == 3
        return 4
end
end

te = Γ.faces[1]
for i = 1:6
    #print(edgee[i],tempor(edgee[i]),abs(CompScienceMeshes.relorientation(edgee[i],te)),"\n")
    @assert tempor(edgee[i]) == abs(CompScienceMeshes.relorientation(edgee[i],te))
end
print("test succesvol")
