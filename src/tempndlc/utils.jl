function relorientation(face::SArray{Tuple{3},Int64,1,3}, tet::SArray{Tuple{4},Int64,1,4})
    v = setdiff(tet,face)
    length(v) == 1 || return 0

    a = something(findfirst(isequal(v[1]), reverse(tet)),0)
    w = sortperm(face)
    b = parity(w)
    return a*(-1)^b
end

function relorientation(edge::SArray{Tuple{2},Int64,1,2}, tet::SArray{Tuple{4},Int64,1,4})
    f = typeof(tet)
    a = zeros(4)
    a[1],a[2],a[3],a[4] = tet[1],tet[2],tet[3],tet[4]
    tet = f([a[4],a[1],a[2],a[3]])
    v = setdiff(tet, edge)
    length(v) == 2 || return 0

    w1 = findfirst(isequal(edge[1]), tet)
    w2 = findfirst(isequal(edge[2]), tet)
    t = tetschoice(w1,w2)
    s = sign(w2-w1)

    return s*t
end

function tetschoice(w1::Int64,w2::Int64)
    a = min(w1,w2)
    b = max(w1,w2)
    if a == 1
        return b-1
    elseif a == 2
        return b+1
    elseif a == 3
        return b+2
    end
end

function relorientation(face::SArray{Tuple{2},Int64,1,2}, simplex::SArray{Tuple{3},Int64,1,3})

    v = setdiff(simplex, face)
    length(v) == 1 || return 0

    # find the position of the missing vertex
    v = v[1]
    #i = Base.findfirst(simplex, v)
    i = something(findfirst(isequal(v), simplex),0)
    s = (-1)^(i-1)

    # remove that vertex from the simplex
    face2 = Array{Int}(undef,length(simplex)-1)
    for j in 1 : i-1
        face2[j] = simplex[j]
    end
    for j in i : length(simplex)-1
        face2[j] = simplex[j+1]
    end

    # get the permutation that maps face to face2
    #p = indexin(face, face2)
    #p = [ findfirst(face2,v) for v in face ]
    p = [something(findfirst(isequal(v),face2),0) for v in face]

    return s * levicivita(p) * i
end
