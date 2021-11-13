using DrWatson
@quickactivate "MIT6.339"

using PyCall, LaTeXStrings
pplt = pyimport("proplot")

include(srcdir("prj1qn2functions.jl"))

m   = MultiGridModel(N=49,methodR="JC",methodP="JC",methodC="JC")
_,ϵ1 = iterateMultiGrid(m,niter=200,ω=0.5,ν1=1,ν2=1,νc=2,nν=2);
_,ϵ2 = iterateMultiGrid(m,niter=200,ω=0.8,ν1=1,ν2=1,νc=2,nν=2);
_,ϵ3 = iterateMultiGrid(m,niter=200,ω=0.5,ν1=1,ν2=1,νc=2,nν=3);
_,ϵ4 = iterateMultiGrid(m,niter=200,ω=0.8,ν1=1,ν2=1,νc=2,nν=3);
_,ϵ5 = iterateMultiGrid(m,niter=200,ω=0.5,ν1=1,ν2=1,νc=2,nν=4);
_,ϵ6 = iterateMultiGrid(m,niter=200,ω=0.8,ν1=1,ν2=1,νc=2,nν=4);

pplt.close(); fig,axs = pplt.subplots(axwidth=2,ncols=3)

axs[1].plot(ϵ1)
axs[1].plot(ϵ2)
axs[2].plot(ϵ3)
axs[2].plot(ϵ4)
axs[3].plot(ϵ5)
axs[3].plot(ϵ6)
# axs[1].contourf(reshape(ϕ,97,97)')

for ax in axs
    ax.format(
        ylim=(10^-15,10^10),yscale="log",yformatter="log",
        ylabel=L"$L_2$",xlabel="Number of Iterations"
    )
end

fig.savefig(plotsdir("prj1_qn2c_testmultigrid.png"),transparent=false,dpi=250)