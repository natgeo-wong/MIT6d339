using DrWatson
@quickactivate "MIT6.339"
using SpecialFunctions

using PyCall, LaTeXStrings
pplt = pyimport("proplot")

# include(srcdir("prj3","quadrature.jl"))
include(srcdir("prj3","creatematrices.jl"))

# upn = upotnumeric()

x = -5:0.01:5; nx = length(x)
y = -2.5:0.01:2.5; ny = length(y)
upa = zeros(nx,ny)
upn = zeros(nx,ny)
lvls = vcat(-10,-5,-2,-1,-0.5,-0.2,-0.1,0.1,0.2,0.5,1,2,5,10)/10

arr = [[0,2,3],[1,4,5],[0,6,7]]
pplt.close(); fig,axs = pplt.subplots(arr,axwidth=2,aspect=2)

A,ψ,xc,xt = ExteriorMatrices2(5); σn = A \ ψ
for iy = 1 : ny, ix = 1 : nx
    if (x[ix]^2 + y[iy]^2) >= 1
        upa[ix,iy] = upotanalyticqn1([x[ix],y[iy]])
        upn[ix,iy] = upotnumeric(σn,xc,[x[ix],y[iy]])
    else
        upa[ix,iy] = NaN
        upn[ix,iy] = NaN
    end
end
axs[2].contourf(x,y,upn',levels=lvls,extend="both")
axs[3].contourf(x,y,upn'.-upa',levels=lvls,extend="both")

A,ψ,xc,xt = ExteriorMatrices2(10); σn = A \ ψ
for iy = 1 : ny, ix = 1 : nx
    if (x[ix]^2 + y[iy]^2) >= 1
        upa[ix,iy] = upotanalyticqn1([x[ix],y[iy]])
        upn[ix,iy] = upotnumeric(σn,xc,[x[ix],y[iy]])
    else
        upa[ix,iy] = NaN
        upn[ix,iy] = NaN
    end
end
axs[4].contourf(x,y,upn',levels=lvls,extend="both")
axs[4].scatter(xc[:,1],xc[:,2],s=3)
axs[4].scatter(xt[:,1],xt[:,2],s=3)
axs[5].contourf(x,y,upn'.-upa',levels=lvls,extend="both")

A,ψ,xc,xt = ExteriorMatrices2(20); σn = A \ ψ
for iy = 1 : ny, ix = 1 : nx
    if (x[ix]^2 + y[iy]^2) >= 1
        upa[ix,iy] = upotanalyticqn1([x[ix],y[iy]])
        upn[ix,iy] = upotnumeric(σn,xc,[x[ix],y[iy]])
    else
        upa[ix,iy] = NaN
        upn[ix,iy] = NaN
    end
end
c = axs[1].contourf(x,y,upa',levels=lvls,extend="both")
axs[6].contourf(x,y,upn',levels=lvls,extend="both")
axs[7].contourf(x,y,upn'.-upa',levels=lvls,extend="both")

axs[1].format(ltitle=L"(a) Analytical $u$",xlabel="x",ylabel="y")
axs[2].format(ltitle=L"(b) Numerical $\hat{u}$",ultitle="n = 5")
axs[3].format(ltitle=L"(c) Difference $\hat{u} - u$")
axs[4].format(ultitle="n = 10")
axs[6].format(ultitle="n = 20",xlabel="x")

fig.colorbar(c,loc="r")
fig.savefig(plotsdir("prj3_qn1d_field.png"),transparent=false,dpi=200)