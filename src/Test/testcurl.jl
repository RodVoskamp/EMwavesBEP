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

X = BEAST.nedelecc3d(Γ)
ref = refspace(X)
ref2 = BEAST.NDLCDRefSpace{Float64}()
t = curl(X)

ptch = chart(Γ, first(cells(Γ)))
ctrd = neighborhood(ptch, [0.3,0.3,0.3])

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

for (ill,k) in enumerate(cells(skeleton(Γ,1)))
    i = tempor(k)
    #print(i)
    #print(ref(ctrd)[i][2], "=")
    a = t.fns[ill][1]
    b = t.fns[ill][2]
    a2 = a.coeff*(ref2(ctrd)[a.refid].value)
    b2 = b.coeff*(ref2(ctrd)[b.refid].value)
    #print(a2+b2)
    @assert (BEAST.norm(ref(ctrd)[i][2]-a2-b2)<0.01)
end
print("test succesvol")
