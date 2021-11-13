using DrWatson
@quickactivate "MIT6.339"
using Printf
using StatsBase

using PyCall, LaTeXStrings
pplt = pyimport("proplot")

include(srcdir("prj1qn1functions.jl"))

Q = zeros(51,50)
MI = zeros(51,50)

for ih = 1 : 50, ib = 1 : 51

    b = (ib-1) * 0.02
    h = ih * 0.02
    m = ChannelModel(b,h,3.,21); _,Q[ib,ih] = SolveChannelModel(m)
    y = h^2 / (b + 2*m.d)
    MI[ib,ih] = y^2 * b * 0.05 + 2 * m.d * 0.05 *(h^2 - 3 * h * y + 3 * y^2) / 3

end

bvec = 0   : 0.02 : 1
hvec = 0.02 : 0.02 : 1

pplt.close(); fig,axs = pplt.subplots(ncols=3,sharey=0,sharex=0,axwidth=1.5)

c = axs[1].contourf(bvec,hvec,Q',levels=(6:60)/400,extend="both")
axs[1].colorbar(c,loc="r")
axs[1].format(ltitle=L"(a) Flow Rate $\hat{Q}$",xlabel="b",ylabel="h")

c = axs[2].contourf(bvec,hvec,MI',levels=(3:30)/2000,extend="both")
axs[2].colorbar(c,loc="r")
axs[2].format(ltitle="(b) Moment of Inertia I",xlabel="b",ylabel="h")

axs[3].scatter(Q[:],MI[:],s=5)
axs[3].format(ltitle="(c) Pareto Optimal Frontier",xlabel=L"Flow Rate $\hat{Q}$",ylabel="Moment of Inertia I")

fig.savefig(plotsdir("prj1_qn1g.png"),transparent=false,dpi=250)