using BEAST
using PyPlot
using CompScienceMeshes
using StaticArrays
using Plots

pyplot()

off = 0
pp = point(0+off,0+off,0)
a = point(1+off,0+off,0)
b = point(0+off,1+off,0)

trip = simplex(a,b,pp)
x = BEAST.NDRefSpace{Float64}()

pygui(true)

fig = figure()
ax = fig.gca()

s = 12
for a in 1:s+1
    for b in 1:s+1
        tu,tv = off+(a-1)/s,off+(b-1)/s
        if tu+tv<=2*off+1
            nbd = neighborhood(trip, [tu,tv])
            q = carttobary(trip,x(nbd)[3][1])
            #q2 = sqrt(q[1]^2+q[2]^2)
            #q1,q2 = q[1]-off,q[2]-off
            ax.quiver(tu,tv, q[1],q[2], scale = 8)
            ax.set_aspect("equal")
        end
    end
end
plt.plot([0+off, 1+off], [0+off, 0+off], color="red")
plt.plot([0+off, 0+off], [0+off, 1+off], color="red")
plt.plot([0+off, 1+off], [1+off, 0+off], color="red")

fig.canvas.draw()
fig.show()
