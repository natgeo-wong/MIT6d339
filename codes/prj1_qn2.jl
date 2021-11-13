using DrWatson
@quickactivate "MIT6.339"

using PyCall, LaTeXStrings
pplt = pyimport("proplot")

include(srcdir("prj1qn2functions.jl"))

m = GaussSeidelModel(N=25)
ϕn,lgn,rgn,ugn,bgn = iterateModel(m,niter=1000,ω=1.5)
ϕr,lgr,rgr,ugr,bgr = directsolveModel(m)
ϕn = reshape(ϕn,m.N,m.N)
ϕr = reshape(ϕr,m.N,m.N)

x = (0 : (m.N-1)) / (m.N-1) * 6

pplt.close(); fig,axs = pplt.subplots(axwidth=2,ncols=3)

c = axs[1].contourf(x,x,ϕn',levels=(0:10)/500,extend="both")
axs[1].colorbar(c,loc="r")

fig.savefig(plotsdir("test.png"),transparent=false,dpi=250)