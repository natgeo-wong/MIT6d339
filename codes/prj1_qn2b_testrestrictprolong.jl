using DrWatson
@quickactivate "MIT6.339"

using PyCall, LaTeXStrings
pplt = pyimport("proplot")

include(srcdir("prj1qn2functions.jl"))

m  = GaussSeidelModel(N=25); ϕ,_,_,_,_,_ = iterateModel(m,niter=750,ω=1.5);
rA,rN = restrictionMatrix(m.N)
pA,pN = prolongateMatrix(m.N)
ϕr = rA * ϕ
ϕp = pA * ϕ

x  = (0 : (m.N-1)) / (m.N-1)
xr = (0 : (m.N-1)/2) / ((m.N-1)/2)
xp = (0 : (m.N-1)*2) / ((m.N-1)*2)

ϕ  = reshape(ϕ,25,25)
ϕr = reshape(ϕr,13,13)
ϕp = reshape(ϕp,49,49)

pplt.close(); fig,axs = pplt.subplots(axwidth=2,ncols=3)

c = axs[1].contourf(x,x,ϕ'*1e3,levels=(2:10)*2.5,extend="both")
axs[2].contourf(xr,xr,ϕr'*1e3,levels=(2:10)*2.5,extend="both")
axs[3].contourf(xp,xp,ϕp'*1e3,levels=(2:10)*2.5,extend="both")

axs[1].format(ltitle="(a) Initial Solution (N=25)")
axs[2].format(ltitle="(b) Restricted Solution (N=13)")
axs[3].format(ltitle="(c) Prolonged Solution (N=49)")

for ax in axs
    ax.plot([1,2,2,1,1]/6,[1,1,2,2,1]/6,c="k")
    ax.plot([2,3,3,2,2]/6,[3,3,4,4,3]/6,c="k")
    ax.plot([4,5,5,4,4]/6,[4,4,5,5,4]/6,c="k")
    ax.plot([4,5,5,4,4]/6,[2,2,3,3,2]/6,c="k")
    ax.format(xlocator=0:1/6:1,xticklabels=[],ylocator=0:1/6:1,yticklabels=[])
end

fig.format(suptitle="Prolongate and Rectriction Tests")
fig.colorbar(c,loc="r",label=L"$\phi$ / 10$^{-3}$")
fig.savefig(plotsdir("prj1_qn2c_testrestrictprolong.png"),transparent=false,dpi=250)