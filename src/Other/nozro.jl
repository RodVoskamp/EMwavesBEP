function nozero(X)
    a,b = size(X)
    c = zeros(a)*im
    d = zeros(b)*im
    k1 = a
    k2 = b
    for i = 1:a
        c[i] = sum(X[i,:])
        if c[i] == 0
            k1 -= 1
        end
    end
    print(k1)
    for i = 1:b
        d[i] = sum(X[:,i])
        if d[i] == 0
            k2 -= 1
        end
    end
    l1,l2 = k1,k2
    X2 = zeros(a,b)*im
    for i = 1:a
        j = a+1-i
        if c[j] != 0
            X2[k1,:] = X[j, :]
            k1 -= 1
        end
    end
    X3 = zeros(a,b)*im
    for i = 1:b
        j = b+1-i
        if d[j] != 0
            X3[:,k2] = X2[:,j]
            k2 -= 1
        end
    end
    return X3[1:l1,1:l2]
end

function nozero(X,u)
    a,b = size(X)
    c = zeros(a)*im
    k1 = a
    for i = 1:a
        c[i] = sum(X[i,:])
        if c[i] == 0
            k1 -= 1
        end
    end
    u2 = zeros(a)*im
    l1 = k1
    for i = 1:a
        j = a+1-i
        if c[j] != 0
            u2[k1] = u[j]
            k1 -= 1
        end
    end
    return u2[1:l1]
end
