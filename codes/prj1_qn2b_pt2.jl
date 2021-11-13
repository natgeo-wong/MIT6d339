using DrWatson
@quickactivate "MIT6.339"

using PyCall, LaTeXStrings
pplt = pyimport("proplot")

include(srcdir("prj1qn2functions.jl"))

bclr = pplt.get_colors("Blues",22)
gclr = pplt.get_colors("Greens",22)

nJC = 10000; mJC = JacobiModel(N=25)
nGS = 2000;  mGS = GaussSeidelModel(N=25)

pplt.close(); fig,axs = pplt.subplots(axwidth=2.5,ncols=2)

for iω = 1 : 20
    _,εJC,_,_,_,_ = iterateModel(mJC,niter=nJC,ω=iω/20)
    axs[1].plot(εJC,c=bclr[iω+1])
end
for iω = 1 : 40
    _,εGS,_,_,_,_ = iterateModel(mGS,niter=nGS,ω=iω/20)
    if iω <= 20
          axs[2].plot(εGS,c=bclr[iω+1],label="$(iω/20)",legend="r",legend_kw=Dict("frame"=>false,"ncol"=>3))
    else; axs[2].plot(εGS,c=gclr[42-iω],label="$(iω/20)",legend="r")
    end
end

axs[1].format(ltitle="(a) Jacobi Iterative Method")
axs[2].format(ltitle="(b) Gauss-Seidel Iterative Method")

for ax in axs
    ax.format(
        ylim=(10^-15,10^10),yscale="log",yformatter="log",
        ylabel=L"$L_2$",xlabel="Number of Iterations"
    )
end

fig.savefig(plotsdir("prj1_qn2b_converge.png"),transparent=false,dpi=250)