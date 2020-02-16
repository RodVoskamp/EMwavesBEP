using BEAST
using CompScienceMeshes
using StaticArrays

tetrs = CompScienceMeshes.tetmeshsphere(1.0, 0.40)
X = BEAST.nedelecc3d(tetrs, skeleton(tetrs,1))
Id = BEAST.Identity()

function f2(p)
    1
end

MP = BEAST.Multiplicative(f2)
Y = curl(X)
A1 = assemble(Id, Y, Y)
A2 = assemble(MP, X, X)
