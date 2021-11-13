using DrWatson
@quickactivate "MIT6.339"
using Printf
using PyCall, LaTeXStrings
pplt = pyimport("proplot")

include(srcdir("prj1qn1functions.jl"))

pplt.close(); fig,axs = pplt.subplots(ncols=3,axwidth=2)

m = ChannelModel(0.5,1.,3.,21)
u,Q = SolveChannelModel(m)
x,y = ξη2xy(m)
axs[1].tricontourf(x[:],y[:],u,cmap="Blues",levels=(0:2:18)/100,extend="both")
axs[1].plot([m.b/2,m.b/2+m.a],[0,m.h],lw=0.5,c="k",linestyle="--")
axs[1].format(xlabel=L"x",ylabel=L"y",ltitle="(a) N = $(m.n)",rtitle="Q = $(@sprintf("%.5f",Q))")

# m = ChannelModel(0.5,1.,3.,51)
# u,Q = SolveChannelModel(m)
# x,y = ξη2xy(m)
# axs[2].tricontourf(x[:],y[:],u,cmap="Blues",levels=(0:2:18)/100,extend="both")
# axs[2].plot([m.b/2,m.b/2+m.a],[0,m.h],lw=0.5,c="k",linestyle="--")
# axs[2].format(xlabel=L"x",ylabel=L"y",ltitle="(a) N = $(m.n)",rtitle="Q = $(@sprintf("%.5f",Q))")

# m = ChannelModel(0.5,1.,3.,101)
# u,Q = SolveChannelModel(m)
# x,y = ξη2xy(m)
# c = axs[3].tricontourf(x[:],y[:],u,cmap="Blues",levels=(0:2:18)/100,extend="both")
# axs[3].plot([m.b/2,m.b/2+m.a],[0,m.h],lw=0.5,c="k",linestyle="--")
# axs[3].colorbar(c,loc="r",label="u")
# axs[3].format(xlabel=L"x",ylabel=L"y",ltitle="(a) N = $(m.n)",rtitle="Q = $(@sprintf("%.5f",Q))")

fig.savefig(plotsdir("prj1_qn1_numsol.png"),transparent=false,dpi=250)