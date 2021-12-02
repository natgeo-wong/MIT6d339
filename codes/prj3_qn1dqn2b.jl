using DrWatson
@quickactivate "MIT6.339"
using SpecialFunctions

using PyCall, LaTeXStrings
pplt = pyimport("proplot")

include(srcdir("prj3","creatematrices.jl"))

nvec = vcat(5:10,12:2:100,120:20:1000,1200:200:5000)
ln = length(nvec)
E1_1d = zeros(ln); E1_2a = zeros(ln); E1_3c = zeros(ln)
E2_1d = zeros(ln); E2_2a = zeros(ln); E2_3c = zeros(ln)
E3_1d = zeros(ln); E3_2a = zeros(ln); E3_3c = zeros(ln)

for ni in 1 : ln
    A,ψ,xc,xt = InteriorMatrices1(nvec[ni]); σn = A \ ψ
    nx = size(xc,1)
    ua = zeros(nx)
    un = zeros(nx)
    
    for ix = 1 : nx
        ua[ix] = upotanalyticqn2(xt[ix,:] * (0.75))
        un[ix] = upotnumeric(σn,xc,xt[ix,:] * (0.75))
    end
    E1_2a[ni] = sum(abs.(un.-ua)) / nvec[ni]
    
    for ix = 1 : nx
        ua[ix] = upotanalyticqn2(xt[ix,:] * (1 - 1/16))
        un[ix] = upotnumeric(σn,xc,xt[ix,:] * (1 - 1/16))
    end
    E2_2a[ni] = sum(abs.(un.-ua)) / nvec[ni]
    
    for ix = 1 : nx
        ua[ix] = upotanalyticqn2(xt[ix,:] * (1 - 1/256))
        un[ix] = upotnumeric(σn,xc,xt[ix,:] * (1 - 1/256))
    end
    E3_2a[ni] = sum(abs.(un.-ua)) / nvec[ni]
end

for ni in 1 : ln
    A,ψ,xc,xt = ExteriorMatrices2(nvec[ni]); σn = A \ ψ
    nx = size(xc,1)
    ua = zeros(nx)
    un = zeros(nx)
    
    for ix = 1 : nx
        ua[ix] = upotanalyticqn1(xt[ix,:] * (1.25))
        un[ix] = upotnumeric(σn,xc,xt[ix,:] * (1.25))
    end
    E1_1d[ni] = sum(abs.(un.-ua)) / nvec[ni]
    
    for ix = 1 : nx
        ua[ix] = upotanalyticqn1(xt[ix,:] * (1 + 1/16))
        un[ix] = upotnumeric(σn,xc,xt[ix,:] * (1 + 1/16))
    end
    E2_1d[ni] = sum(abs.(un.-ua)) / nvec[ni]
    
    for ix = 1 : nx
        ua[ix] = upotanalyticqn1(xt[ix,:] * (1 + 1/256))
        un[ix] = upotnumeric(σn,xc,xt[ix,:] * (1 + 1/256))
    end
    E3_1d[ni] = sum(abs.(un.-ua)) / nvec[ni]
end

a = 1

for ni in 1 : ln
    A,ψ,xc,xt = ExteriorMatrices3(nvec[ni],a=a); σn = A \ ψ
    nx = size(xc,1)
    ua = zeros(nx)
    un = zeros(nx)
    
    for ix = 1 : nx
        ua[ix] = upotanalyticqn1([xt[ix,1],xt[ix,2]*a] * (1 + 1/4))
        un[ix] = upotnumericθ(σn,xc,[xt[ix,1],xt[ix,2]*a] * (1 + 1/4),a=a)
    end
    E1_3c[ni] = sum(abs.(un.-ua)) / nvec[ni]
    
    for ix = 1 : nx
        ua[ix] = upotanalyticqn1([xt[ix,1],xt[ix,2]*a] * (1 + 1/16))
        un[ix] = upotnumericθ(σn,xc,[xt[ix,1],xt[ix,2]*a] * (1 + 1/16),a=a)
    end
    E2_3c[ni] = sum(abs.(un.-ua)) / nvec[ni]
    
    for ix = 1 : nx
        ua[ix] = upotanalyticqn1([xt[ix,1],xt[ix,2]*a] * (1 + 1/256))
        un[ix] = upotnumericθ(σn,xc,[xt[ix,1],xt[ix,2]*a] * (1 + 1/256),a=a)
    end
    E3_3c[ni] = sum(abs.(un.-ua)) / nvec[ni]
end

pplt.close(); fig,axs = pplt.subplots(ncols=3,axwidth=2)

lgd = Dict("frame"=>false,"ncol"=>1)

axs[1].plot(nvec,E1_1d)
axs[1].plot(nvec,E2_1d)
axs[1].plot(nvec,E3_1d)
axs[1].plot(10. .^[0,4],10. .^[0,-8],lw=0.5,c="k",linestyle=":")

axs[2].plot(nvec,E1_2a)
axs[2].plot(nvec,E2_2a)
axs[2].plot(nvec,E3_2a)
axs[2].plot(10. .^[0,4],10. .^[0,-8],lw=0.5,c="k",linestyle=":")

axs[3].plot(nvec,E1_3c,label=L"$\delta$r = $\pm\frac{1}{4}$",legend="r",legend_kw=lgd)
axs[3].plot(nvec,E2_3c,label=L"$\delta$r = $\pm\frac{1}{16}$",legend="r")
axs[3].plot(nvec,E3_3c,label=L"$\delta$r = $\pm\frac{1}{256}$",legend="r")
axs[3].plot(10. .^[0,4],10. .^[0,-8],label=L"$\alpha$ = 2",legend="r",lw=0.5,c="k",linestyle=":")
axs[3].plot(10. .^[0,4],10. .^[0,12],label=L"$\alpha$ = 3",legend="r",lw=0.5,c="k",linestyle="--")

axs[1].format(ltitle="(a) Exterior Neumann (Q1d)")
axs[2].format(ltitle="(b) Interior Neumann (Q2b)")
axs[3].format(ltitle="(c) Ellipse Problem (Q3c)")

for ax in axs
    ax.format(xscale="log",yscale="log",xformatter="log",yformatter="log")
    ax.format(xlim=(5,5000),ylim=10. .^(-8,0),xlabel="n",ylabel=L"|u-\hat{u}|")
end

fig.savefig(plotsdir("prj3_converge.png"),transparent=false,dpi=200)