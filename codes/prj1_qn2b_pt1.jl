using DrWatson
@quickactivate "MIT6.339"

using PyCall, LaTeXStrings
pplt = pyimport("proplot")

include(srcdir("prj1qn2functions.jl"))

mJC = JacobiModel(N=25)
ϕJ_1,_,_,_,_,_ = iterateModel(mJC,niter=200,ω=1.5); eJ_1 = mJC.A * ϕJ_1 .- mJC.f
ϕJ_2,_,_,_,_,_ = iterateModel(mJC,niter=200,ω=1.0); eJ_2 = mJC.A * ϕJ_2 .- mJC.f
ϕJ_3,_,_,_,_,_ = iterateModel(mJC,niter=200,ω=0.3); eJ_3 = mJC.A * ϕJ_3 .- mJC.f

mGS = GaussSeidelModel(N=25)
ϕGS1,_,_,_,_,_ = iterateModel(mGS,niter=200,ω=1.5); eGS1 = mGS.A * ϕGS1 .- mGS.f
ϕGS2,_,_,_,_,_ = iterateModel(mGS,niter=200,ω=1.0); eGS2 = mGS.A * ϕGS2 .- mGS.f
ϕGS3,_,_,_,_,_ = iterateModel(mGS,niter=200,ω=0.3); eGS3 = mGS.A * ϕGS3 .- mGS.f

ϕr = mJC.A \ Vector(mJC.f)

x = (0 : (mJC.N-1)) / (mJC.N-1)

arr = [[0,2,3,4],[1,2,3,4],[1,5,6,7],[0,5,6,7]]
pplt.close(); fig,axs = pplt.subplots(arr,axwidth=1.5,hspace=0,wspace=[3,0,0])

c = axs[1].contourf(x,x,reshape(ϕr,mJC.N,mJC.N)'*1e3,levels=(2:10)*2.5,extend="both")

axs[2].contourf(x,x,reshape(ϕJ_1,mJC.N,mJC.N)'*1e3,levels=(2:10)*2.5,extend="both")
axs[3].contourf(x,x,reshape(ϕJ_2,mJC.N,mJC.N)'*1e3,levels=(2:10)*2.5,extend="both")
axs[4].contourf(x,x,reshape(ϕJ_3,mJC.N,mJC.N)'*1e3,levels=(2:10)*2.5,extend="both")

axs[5].contourf(x,x,reshape(ϕGS1,mGS.N,mGS.N)'*1e3,levels=(2:10)*2.5,extend="both")
axs[6].contourf(x,x,reshape(ϕGS2,mGS.N,mGS.N)'*1e3,levels=(2:10)*2.5,extend="both")
axs[7].contourf(x,x,reshape(ϕGS3,mGS.N,mGS.N)'*1e3,levels=(2:10)*2.5,extend="both")

axs[1].format(ultitle="(a) Direct Solve")
axs[2].format(ultitle=L"(b) Jacobi, $\omega$ = 1.5")
axs[3].format(ultitle=L"(c) Jacobi, $\omega$ = 1.0")
axs[4].format(ultitle=L"(d) Jacobi, $\omega$ = 0.3")
axs[5].format(ultitle=L"(e) G-S, $\omega$ = 1.5")
axs[6].format(ultitle=L"(f) G-S, $\omega$ = 1.0")
axs[7].format(ultitle=L"(g) G-S, $\omega$ = 0.3")

for ax in axs
    ax.plot([1,2,2,1,1]/6,[1,1,2,2,1]/6,c="k")
    ax.plot([2,3,3,2,2]/6,[3,3,4,4,3]/6,c="k")
    ax.plot([4,5,5,4,4]/6,[4,4,5,5,4]/6,c="k")
    ax.plot([4,5,5,4,4]/6,[2,2,3,3,2]/6,c="k")
    ax.format(xlocator=0:1/6:1,xticklabels=[],ylocator=0:1/6:1,yticklabels=[])
end

fig.format(suptitle="N = 25 | 200 Iterations | Source Blocks (1,7,14,16)")
fig.colorbar(c,loc="r",label=L"$\phi$ / 10$^{-3}$")
fig.savefig(plotsdir("prj1_qn2_solutions.png"),transparent=false,dpi=250)