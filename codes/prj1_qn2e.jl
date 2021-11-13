using DrWatson
@quickactivate "MIT6.339"

using PyCall, LaTeXStrings
pplt = pyimport("proplot")

include(srcdir("prj1qn2functions.jl"))

m = MultiGridModel(N=49)

_,ϵ1_02 = iterateMultiGrid(m,niter=250,ω=0.5,ν1=1,ν2=1,νc=2,nν=2)
_,ϵ1_03 = iterateMultiGrid(m,niter=250,ω=0.8,ν1=1,ν2=1,νc=2,nν=3)
_,ϵ1_04 = iterateMultiGrid(m,niter=250,ω=0.5,ν1=1,ν2=1,νc=2,nν=4)
_,ϵ2_02 = iterateMultiGrid(m,niter=250,ω=0.5,ν1=1,ν2=1,νc=5,nν=2)
_,ϵ2_03 = iterateMultiGrid(m,niter=250,ω=0.8,ν1=1,ν2=1,νc=5,nν=3)
_,ϵ2_04 = iterateMultiGrid(m,niter=250,ω=0.5,ν1=1,ν2=1,νc=5,nν=4)
_,ϵ3_02 = iterateMultiGrid(m,niter=250,ω=0.5,ν1=1,ν2=1,νc=10,nν=2)
_,ϵ3_03 = iterateMultiGrid(m,niter=250,ω=0.8,ν1=1,ν2=1,νc=10,nν=3)
_,ϵ3_04 = iterateMultiGrid(m,niter=250,ω=0.5,ν1=1,ν2=1,νc=10,nν=4)

_,ϵ4_02 = iterateMultiGrid(m,niter=250,ω=0.5,ν1=2,ν2=2,νc=2,nν=2)
_,ϵ4_03 = iterateMultiGrid(m,niter=250,ω=0.8,ν1=2,ν2=2,νc=2,nν=3)
_,ϵ4_04 = iterateMultiGrid(m,niter=250,ω=0.5,ν1=2,ν2=2,νc=2,nν=4)
_,ϵ5_02 = iterateMultiGrid(m,niter=250,ω=0.5,ν1=2,ν2=2,νc=5,nν=2)
_,ϵ5_03 = iterateMultiGrid(m,niter=250,ω=0.8,ν1=2,ν2=2,νc=5,nν=3)
_,ϵ5_04 = iterateMultiGrid(m,niter=250,ω=0.5,ν1=2,ν2=2,νc=5,nν=4)
_,ϵ6_02 = iterateMultiGrid(m,niter=250,ω=0.5,ν1=2,ν2=2,νc=10,nν=2)
_,ϵ6_03 = iterateMultiGrid(m,niter=250,ω=0.8,ν1=2,ν2=2,νc=10,nν=3)
_,ϵ6_04 = iterateMultiGrid(m,niter=250,ω=0.5,ν1=2,ν2=2,νc=10,nν=4)

pplt.close(); fig,axs = pplt.subplots(axwidth=2,ncols=3,nrows=2)

lgd = Dict("ncols"=>3,"frame"=>false)

axs[1].plot(ϵ1_02); axs[2].plot(ϵ2_02); axs[3].plot(ϵ2_02)
axs[1].plot(ϵ1_03); axs[2].plot(ϵ2_03); axs[3].plot(ϵ2_03)
axs[1].plot(ϵ1_04); axs[2].plot(ϵ2_04); axs[3].plot(ϵ2_04)

axs[4].plot(ϵ4_02); axs[5].plot(ϵ5_02,label=L"$n\nu$ = 2",legend="b",legend_kw=lgd)
axs[4].plot(ϵ4_03); axs[5].plot(ϵ5_03,label=L"$n\nu$ = 3",legend="b")
axs[4].plot(ϵ4_04); axs[5].plot(ϵ5_04,label=L"$n\nu$ = 4",legend="b")

axs[6].plot(ϵ6_02)
axs[6].plot(ϵ6_03)
axs[6].plot(ϵ6_04)

axs[1].format(ultitle=L"(a) $\nu_1=\nu_2=1$, $\nu_c=2$")
axs[2].format(ultitle=L"(b) $\nu_1=\nu_2=1$, $\nu_c=5$")
axs[3].format(ultitle=L"(c) $\nu_1=\nu_2=1$, $\nu_c=10$")
axs[4].format(ultitle=L"(a) $\nu_1=\nu_2=2$, $\nu_c=2$")
axs[5].format(ultitle=L"(b) $\nu_1=\nu_2=2$, $\nu_c=5$")
axs[6].format(ultitle=L"(c) $\nu_1=\nu_2=2$, $\nu_c=10$")

for ax in axs
    ax.format(
        ylim=(10^-11,0.1),yscale="log",yformatter="log",
        ylabel=L"$L_2$",xlabel="Number of Cycles",
        suptitle="2-, 3-, and 4-Grid Refinement V-Cycles"
    )
end

fig.savefig(plotsdir("prj1_qn2e.png"),transparent=false,dpi=250)