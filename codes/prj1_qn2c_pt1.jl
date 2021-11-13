using DrWatson
@quickactivate "MIT6.339"

using PyCall, LaTeXStrings
pplt = pyimport("proplot")

include(srcdir("prj1qn2functions.jl"))

mGGG = MultiGridModel()
mGGJ = MultiGridModel(C="JC")
mGJJ = MultiGridModel(P="JC",C="JC")
mGJG = MultiGridModel(P="JC")
mJGG = MultiGridModel(R="JC",)
mJGJ = MultiGridModel(R="JC",C="JC")
mJJG = MultiGridModel(R="JC",P="JC")
mJJJ = MultiGridModel(R="JC",P="JC",C="JC")

_,ϵGGG1 = iterateMultiGrid(mGGG,niter=300,ω=0.5,ν1=1,ν2=1,νc=2,nν=2);
_,ϵGGJ1 = iterateMultiGrid(mGGJ,niter=300,ω=0.5,ν1=1,ν2=1,νc=2,nν=2);
_,ϵGJJ1 = iterateMultiGrid(mGJJ,niter=300,ω=0.5,ν1=1,ν2=1,νc=2,nν=2);
_,ϵGJG1 = iterateMultiGrid(mGJG,niter=300,ω=0.5,ν1=1,ν2=1,νc=2,nν=2);
_,ϵJGG1 = iterateMultiGrid(mJGG,niter=300,ω=0.5,ν1=1,ν2=1,νc=2,nν=2);
_,ϵJGJ1 = iterateMultiGrid(mJGJ,niter=300,ω=0.5,ν1=1,ν2=1,νc=2,nν=2);
_,ϵJJG1 = iterateMultiGrid(mJJG,niter=300,ω=0.5,ν1=1,ν2=1,νc=2,nν=2);
_,ϵJJJ1 = iterateMultiGrid(mJJJ,niter=300,ω=0.5,ν1=1,ν2=1,νc=2,nν=2);

_,ϵGGG2 = iterateMultiGrid(mGGG,niter=300,ω=0.8,ν1=1,ν2=1,νc=2,nν=2);
_,ϵGGJ2 = iterateMultiGrid(mGGJ,niter=300,ω=0.8,ν1=1,ν2=1,νc=2,nν=2);
_,ϵGJJ2 = iterateMultiGrid(mGJJ,niter=300,ω=0.8,ν1=1,ν2=1,νc=2,nν=2);
_,ϵGJG2 = iterateMultiGrid(mGJG,niter=300,ω=0.8,ν1=1,ν2=1,νc=2,nν=2);
_,ϵJGG2 = iterateMultiGrid(mJGG,niter=300,ω=0.8,ν1=1,ν2=1,νc=2,nν=2);
_,ϵJGJ2 = iterateMultiGrid(mJGJ,niter=300,ω=0.8,ν1=1,ν2=1,νc=2,nν=2);
_,ϵJJG2 = iterateMultiGrid(mJJG,niter=300,ω=0.8,ν1=1,ν2=1,νc=2,nν=2);
_,ϵJJJ2 = iterateMultiGrid(mJJJ,niter=300,ω=0.8,ν1=1,ν2=1,νc=2,nν=2);

pplt.close(); fig,axs = pplt.subplots(axwidth=1.5,ncols=4,nrows=2)

lgd = Dict("ncols"=>1,"frame"=>false)

axs[1].plot(ϵGGG1); axs[1].plot(ϵGGG2)
axs[2].plot(ϵGGJ1); axs[2].plot(ϵGGJ2)
axs[3].plot(ϵGJJ1); axs[3].plot(ϵGJJ2)
axs[4].plot(ϵGJG1); axs[4].plot(ϵGJG2)
axs[5].plot(ϵJGG1); axs[5].plot(ϵJGG2)
axs[6].plot(ϵJGJ1); axs[6].plot(ϵJGJ2)
axs[7].plot(ϵJJG1); axs[7].plot(ϵJJG2)
axs[8].plot(ϵJJJ1); axs[8].plot(ϵJJJ2)

axs[1].format(ultitle="(a) GGG")
axs[2].format(ultitle="(b) GGJ")
axs[3].format(ultitle="(c) GJJ")
axs[4].format(ultitle="(d) GJG")
axs[5].format(ultitle="(e) JGG")
axs[6].format(ultitle="(f) JGJ")
axs[7].format(ultitle="(g) JJG")
axs[8].format(ultitle="(h) JJJ")

for ax in axs
    ax.format(
        ylim=(10^-11,0.1),yscale="log",yformatter="log",
        ylabel=L"$L_2$",xlabel="Number of Iterations",
        suptitle="Iterative Method Combinations"
    )
end

fig.savefig(plotsdir("prj1_qn2c_pt1.png"),transparent=false,dpi=250)