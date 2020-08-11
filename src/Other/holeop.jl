t = range(-2,2,length=200)
t2 = range(-2,2,length=200)
pts = [point(s,s2) for s in t, s2 in t2]

function holeop(R)
    A = vec(ones(size(pts)))
    i = 1
    for pt in pts
        if R<sqrt(pt[2]^2+pt[1]^2)< 1
            A[i] = (sqrt(pt[2]^2+pt[1]^2)-R)/(1-R)
        end
        i += 1
    end
    return A
end

AA = zeros(0)
for r in range(0,stop=0.4,length=100)
    global AA
    append!(AA,holeop(r))
end

anim = @animate for i ∈ 1:100
    display(plot(surface(t,t2,AA[length(A)*(i-1)+1:length(A)*i],aspect_ratio=:equal,clim=(0,1),camera=(0,90),size=(500,400),c=ColorGradient([:white,:yellow,:green,:orange,:red,:blue,:grey,:black]))))
end
gif(anim, "holeop.gif", fps = 15)
