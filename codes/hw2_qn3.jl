using DrWatson
@quickactivate "MIT6.339"
using LinearAlgebra
using StatsBase
using PyCall, LaTeXStrings
pplt = pyimport("proplot")

include(srcdir("hw2functions.jl"))

nxini = 5
nitr = 8

clr = pplt.Colors("Blues",nitr+1)

dx = zeros(nitr+1)
p1_1 = zeros(nitr+1); p1_2 = zeros(nitr+1); p1_3 = zeros(nitr+1)
p2_1 = zeros(nitr+1); p2_2 = zeros(nitr+1); p2_3 = zeros(nitr+1)
p∞_1 = zeros(nitr+1); p∞_2 = zeros(nitr+1); p∞_3 = zeros(nitr+1)

pplt.close(); fig,axs = pplt.subplots(ncols=3,sharey=1,axwidth=2)

for j = 0 : nitr

    nx = nxini * 2^(j+1)
    X = (1:(nx-1)) / nx
    A = Tridiagonal(ones(nx-2),ones(nx-1)*-2,ones(nx-2)) * nx^2
    dx[j+1] = 1 / nx

    ua1 = u1.(X); un1 = - A \ f1.((1:(nx-1))/nx)
    ua2 = u2.(X); un2 = - A \ f2.((1:(nx-1))/nx)
    ua3 = u3.(X); un3 = - A \ f3.((1:(nx-1))/nx)

    p1_1[j+1] = p1norm(ua1,un1,nx); p2_1[j+1] = p2norm(ua1,un1); p∞_1[j+1] = p∞norm(ua1,un1)
    p1_2[j+1] = p1norm(ua2,un2,nx); p2_2[j+1] = p2norm(ua2,un2); p∞_2[j+1] = p∞norm(ua2,un2)
    p1_3[j+1] = p1norm(ua3,un3,nx); p2_3[j+1] = p2norm(ua3,un3); p∞_3[j+1] = p∞norm(ua3,un3)

    axs[1].scatter(X,un1,s=5,c=clr[j+1])
    axs[2].scatter(X,un2,s=5,c=clr[j+1])
    axs[3].scatter(X,un3,s=5,c=clr[j+1],legend="r",label="nx = $(nxini * 2^(j+1))",legend_kw=Dict("ncols"=>1,"frame"=>false))

    if j == nitr
        axs[1].plot(X,ua1,lw=0.5,c="k")
        axs[2].plot(X,ua2,lw=0.5,c="k")
        axs[3].plot(X,ua3,lw=0.5,c="k")
    end

    axs[1].format(ltitle="(a) f(x) = (6x+2x^2)e^x")
    axs[2].format(ltitle="(b) f(x) = 2 or f(x) = 0")
    axs[3].format(ltitle="(c) f(x) = 3x")

end

for ax in axs
    ax.format(xlabel=L"x",ylabel="f(x)")
end

fig.savefig(plotsdir("hw2_qn3_numericalsol.png"),dpi=250,transparent=false)

pplt.close(); fig,axs = pplt.subplots(ncols=3,sharey=1,axwidth=2)

axs[1].plot(log.(dx),log.(p1_1))
axs[1].plot(log.(dx),log.(p2_1))
axs[1].plot(log.(dx),log.(p∞_1))
axs[1].format(ltitle="(a) f(x) = (6x+2x^2)e^x")

axs[2].plot(log.(dx),log.(p1_2))
axs[2].plot(log.(dx),log.(p2_2))
axs[2].plot(log.(dx),log.(p∞_2))
axs[2].format(ltitle="(b) f(x) = 2 or f(x) = 0")

axs[3].plot(log.(dx),log.(p1_3),legend="r",label="p = 1",legend_kw=Dict("ncols"=>1,"frame"=>false))
axs[3].plot(log.(dx),log.(p2_3),legend="r",label="p = 2")
axs[3].plot(log.(dx),log.(p∞_3),legend="r",label=L"p = $\infty$")
axs[3].format(ltitle="(c) f(x) = 3x")

for ax in axs
    ax.format(xlim=(-9,-1),xlabel=L"log(\Delta x)",ylabel=L"log || v ||_p^*")
end

fig.savefig(plotsdir("hw2_qn3_pnorm.png"),dpi=250,transparent=false)
