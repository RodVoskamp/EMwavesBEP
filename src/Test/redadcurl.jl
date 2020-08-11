using BEAST
using EMwavesBEP
using CompScienceMeshes
using StaticArrays

v = SArray{Tuple{3},Float64,1,3}
f = SArray{Tuple{4},Int,1,4}
f2 = SArray{Tuple{2},Int,1,2}

o = v([0,0,0])
a = v([1,0,0])
b = v([0,2,0])
c = v([0,0,4])
face = zeros(f,1)
face[1] = f([1,2,3,4])
Γ = Mesh([a,b,c,o], face)

edm = zeros(f2,3)
edm[1] = f2([1,2])
edm[2] = f2([1,3])
edm[3] = f2([1,4])
useedg = Mesh([a,b,c,o], edm)
X = BEAST.nedelecc3d(Γ,useedg)

Y = curl(X)
E, ad = assemblydata(Y)
k = 0
for i = 1:4
    l = 0
for (m,a) in ad[1,i] #println(m," ", a)
    global k += 1
    if i != 1
        l += 1
    end
end
if i != 1
    @assert l == 2
end
#print("\n")
end
@assert k == 6

X = BEAST.nedelecc3d(Γ)
Y = curl(X)
E, ad = assemblydata(Y)
k = 0
for i = 1:4
    l = 0
for (m,a) in ad[1,i] #println(m," ", a)
    global k += 1
    l += 1
end
@assert l == 3
#print("\n")
end

@assert k == 12

print("test succesvol")
