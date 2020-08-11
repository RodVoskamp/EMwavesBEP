using BEAST
using CompScienceMeshes
using StaticArrays
using EMwavesBEP
using Plots
using DelimitedFiles
using LinearAlgebra
using SparseArrays
using IterativeSolvers

Ω = CompScienceMeshes.tetmeshsphere(1,0.2)
X = BEAST.nedelecc3d(Ω)
ttrX = BEAST.ttrace(X,boundary(Ω))
Id = BEAST.Identity()
A22 = assemble(Id,X,X)
A2 = Array(A22)
e = Maxwell3D.planewave(direction=ẑ, polarization=x̂, wavenumber=1.0)
b1 = assemble((n × e) × n,ttrX)

A = [4.0 1; 1 3]
A = real(A2)
print("condition number:", cond(A), "\n")
@assert isposdef(A) == true
b = [1.0, 2]
b = real(b1)

A\b
print(typeof(A),"\n")
@time begin
    c = A\b
end
As = sparse(A)
bs = sparse(b)
cs = zeros(length(bs))
gmres!(cs, As, bs)
print(typeof(As),"\n")
@time begin
    cs = gmres!(cs, As, bs)
end
@assert norm(c-cs)<0.01
