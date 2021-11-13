using DrWatson
@quickactivate "MIT6.339"

using PyCall, LaTeXStrings
pplt = pyimport("proplot")

pplt.close(); f,axs = pplt.subplots(ncols=3,axwidth=2)

axs[1].plot([0,1,0,-1,0],[1,0,-1,0,1],c="k")
axs[1].fill([0,1,0,-1,0],[1,0,-1,0,1])
axs[1].format(ultitle="(a) p = 1")

axs[2].plot(sin.(range(0,stop=2*π,length=201)),cos.(range(0,stop=2*π,length=201)),c="k")
axs[2].fill(sin.(range(0,stop=2*π,length=201)),cos.(range(0,stop=2*π,length=201)))
axs[2].format(ultitle="(b) p = 2")

axs[3].plot([-1,-1,1,1,-1],[-1,1,1,-1,-1],c="k")
axs[3].fill([-1,-1,1,1,-1],[-1,1,1,-1,-1])
axs[3].format(ultitle=L"(c) p = $\infty$")

for ax in axs
    ax.format(xlim=(-2,2),ylim=(-2,2),xlabel="x",ylabel="y")
end

f.savefig(plotsdir("hw2_qn1.png"),dpi=200,transparent=false)