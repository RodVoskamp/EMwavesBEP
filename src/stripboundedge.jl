using BEAST
using CompScienceMeshes

function stripboundedge(m::Mesh{3,4,Float64})
    f = skeleton(Γ,2)
    C = connectivity(f,Γ)

    totfaces = length(f.faces)
    v = SArray{Tuple{3},Float64,1,3}
    A = zeros(v,2*totfaces-4*length(Γ.faces))
    i = 1
    k = 1
    while i < totfaces
        if length(BEAST.nzrange(C,i)) ≈ 1
            A[k] = f.faces[i]
            k += 1
        end
        i += 1
    end

    edges = skeleton(Γ,1)

    Γ2 = Mesh(Γ.vertices,A)
    rem = skeleton(Γ2,1)
    savedges = setdiff(edges.faces,rem.faces)
    return savedges, rem.faces
end
