using DrWatson
@quickactivate "MIT6.339"

using PyCall, LaTeXStrings
pplt = pyimport("proplot")

include(srcdir("prj1qn2functions.jl"))

m = MultiGridModel()

_,ϵ1_02 = iterateMultiGrid(m,niter=250,ω=0.5,ν1=1,ν2=1,νc=2,nν=2);
_,ϵ2_02 = iterateMultiGrid(m,niter=250,ω=0.8,ν1=1,ν2=1,νc=2,nν=2);
_,ϵ1_04 = iterateMultiGrid(m,niter=250,ω=0.5,ν1=1,ν2=1,νc=4,nν=2);
_,ϵ2_04 = iterateMultiGrid(m,niter=250,ω=0.8,ν1=1,ν2=1,νc=4,nν=2);
_,ϵ1_10 = iterateMultiGrid(m,niter=250,ω=0.5,ν1=1,ν2=1,νc=10,nν=2);
_,ϵ2_10 = iterateMultiGrid(m,niter=250,ω=0.8,ν1=1,ν2=1,νc=10,nν=2);
_,ϵ1_20 = iterateMultiGrid(m,niter=250,ω=0.5,ν1=1,ν2=1,νc=20,nν=2);
_,ϵ2_20 = iterateMultiGrid(m,niter=250,ω=0.8,ν1=1,ν2=1,νc=20,nν=2);
_,ϵ1_50 = iterateMultiGrid(m,niter=250,ω=0.5,ν1=1,ν2=1,νc=50,nν=2);
_,ϵ2_50 = iterateMultiGrid(m,niter=250,ω=0.8,ν1=1,ν2=1,νc=50,nν=2);
_,ϵ1_∞  = iterateMultiGrid(m,niter=250,ω=0.5,ν1=1,ν2=1,νc=0,nν=2);
_,ϵ2_∞  = iterateMultiGrid(m,niter=250,ω=0.8,ν1=1,ν2=1,νc=0,nν=2);

pplt.close(); fig,axs = pplt.subplots(axwidth=2,ncols=2)

lgd = Dict("ncols"=>1,"frame"=>false)

axs[1].plot(ϵ1_02); axs[2].plot(ϵ2_02,label=L"$\nu_c$ = 2",legend="r",legend_kw=lgd)
axs[1].plot(ϵ1_04); axs[2].plot(ϵ2_04,label=L"$\nu_c$ = 4",legend="r")
axs[1].plot(ϵ1_10); axs[2].plot(ϵ2_10,label=L"$\nu_c$ = 10",legend="r")
axs[1].plot(ϵ1_20); axs[2].plot(ϵ2_20,label=L"$\nu_c$ = 20",legend="r")
axs[1].plot(ϵ1_50); axs[2].plot(ϵ2_50,label=L"$\nu_c$ = 50",legend="r");
axs[1].plot(ϵ1_∞);  axs[2].plot(ϵ2_∞,label=L"$\nu_c$ = $\infty$",legend="r");

axs[1].format(ultitle=L"(a) $\omega=0.5$")
axs[2].format(ultitle=L"(b) $\omega=0.8$")

for ax in axs
    ax.format(
        ylim=(10^-11,0.1),yscale="log",yformatter="log",
        ylabel=L"$L_2$",xlabel="Number of Iterations",
        suptitle=L"Varying $\nu_c$"
    )
end

fig.savefig(plotsdir("prj1_qn2c_pt2.png"),transparent=false,dpi=250)