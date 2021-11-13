using DrWatson
@quickactivate "MIT6.339"
using LinearAlgebra
using StatsBase
using PyCall, LaTeXStrings
pplt = pyimport("proplot")

include(srcdir("hw2functions.jl"))

nxini = 5
nitr = 8

clr = pplt.Colors("Blues",nitr+2)

dx = zeros(nitr+1)
p2_1 = zeros(nitr+1)
p2_2 = zeros(nitr+1)

pplt.close(); fig,axs = pplt.subplots(ncols=3,sharey=0,sharex=0,axwidth=2)

for j = 0 : nitr

    nx = nxini * 2^(j+1)
    X = (1:nx) / nx
    dx[j+1] = 1 / nx

    A1 = Tridiagonal(vcat(ones(nx-2),-1),vcat(ones(nx-1)*-2,1),ones(nx-1)) * nx^2
    A1[end,:] = A1[end,:] / nx

    A2 = Tridiagonal(vcat(ones(nx-1),0),vcat(ones(nx)*-2,1),ones(nx)) * nx^2
    A2[end,:] = A2[end,:] / nx
    A2 = Array(A2)
    A2[end,end-2] = -10

    ua1 = u1.(X)
    un1 = - A1 \ vcat(f1.((1:(nx-1))/nx),2*exp(1))
    un2 = - A2 \ vcat(f1.((1:nx)/nx),2*exp(1))

    p2_1[j+1] = p2norm(ua1,un1);
    p2_2[j+1] = p2norm(ua1,un2[1:(end-1)]);

    axs[1].scatter(X,un1,s=5,c=clr[j+2])
    axs[2].scatter(X,un2[1:(end-1)],s=5,c=clr[j+2],legend="r",label="nx = $(nxini * 2^(j+1))",legend_kw=Dict("ncols"=>1,"frame"=>false))

    if j == nitr
        axs[1].plot(X,ua1,lw=1,c="k")
        axs[2].plot(X,ua1,lw=1,c="k")
    end

    # axs[1].format(ltitle="(a) f(x) = (6x+2x^2)e^x")
end

axs[3].plot(log.(dx),log.(p2_1),legend="lr",label="1st-order backward",legend_kw=Dict("ncols"=>1,"frame"=>false))
axs[3].plot(log.(dx),log.(p2_2),legend="lr",label="2nd-order central")
axs[3].format(xlim=(-8,-2),xlabel=L"log(\Delta x)",ylabel=L"log || v ||_2^*",ltitle="(c) Order of Convergence")

axs[1].format(xlabel=L"x",ltitle="(a) 1st-order Backward",ylabel="f(x)")
axs[2].format(xlabel=L"x",ltitle="(b) 2nd-order Central",)

fig.savefig(plotsdir("hw2_qn4_numericalsol.png"),dpi=250,transparent=false)