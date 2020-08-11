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

facee = skeleton(Γ,2).faces
te = Γ.faces[1]
for i = 1:4
    #print(setdiff(te,facee[i]),CompScienceMeshes.relorientation(facee[i],te))
    @assert setdiff(te,facee[i])[1] == abs(CompScienceMeshes.relorientation(facee[i],te))
    #print(setdiff(te,facee[i])[1] == abs(CompScienceMeshes.relorientation(facee[i],te)),"\n")
end
"""
#test the curl

X = BEAST.nedelecc3d(Γ)
ref = refspace(X)
ref2 = BEAST.NDLCDRefSpace{Float64}()
t = curl(X)

ptch = chart(Γ, first(cells(Γ)))
ctrd = neighborhood(ptch, [0.3,0.3,0.3])

for (i,k) in enumerate(cells(skeleton(Γ,1)))
    #print(ref(ctrd)[i][2], "=")
    a = t.fns[i][1]
    b = t.fns[i][2]
    a2 = a.coeff*(ref2(ctrd)[a.refid].value)
    b2 = b.coeff*(ref2(ctrd)[b.refid].value)
    #print(a2+b2)
    @assert BEAST.norm(ref(ctrd)[i][2]-a2-b2)<0.01
    #print(BEAST.norm(ref(ctrd)[i][2]-a2-b2)<0.01,"\n")
end


#ntr = BEAST.ntrace(t,skeleton(Γ,2))
#ttr = BEAST.ttrace(t,skeleton(Γ,2))
"""

print("test succesvol")
