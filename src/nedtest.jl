using BEAST
using PyPlot
using CompScienceMeshes
using StaticArrays
using Plots

pyplot()

d = dirname(pathof(BEAST))
include(joinpath(d,"../src/bases/local/ndlcdlocal.jl"))

o = point(0,0,0)
x = point(1,0,0)
y = point(0,1,0)
z = point(0,0,1)

tet = simplex(o,x,y,z)
x = NDLCDRefSpace{Float64}()

pygui(true)

fig = figure()
ax = fig.gca(projection="3d")

for a in 1:11
    for b in 1:11
        for c in 1:11
            tu,tv,tw = (a-1.0)/10,(b-1.0)/10,(c-1.0)/10
            if tu+tv+tw<=1
                nbd = neighborhood(tet, [tu,tv,tw])
                q = x(nbd)[2][1]/5
                ax.quiver(tu,tv,tw, q[1],q[2],q[3])
            end
        end
    end
end
fig.canvas.draw()
