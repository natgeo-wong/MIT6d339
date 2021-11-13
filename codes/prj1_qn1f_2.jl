using DrWatson
@quickactivate "MIT6.339"
using Printf
using StatsBase

using PyCall, LaTeXStrings
pplt = pyimport("proplot")

include(srcdir("prj1qn1functions.jl"))

pplt.close(); fig,axs = pplt.subplots(nrows=2,ncols=3,axwidth=1.5,sharex=0,sharey=0)

m1 = ChannelModel(1.,1.,3.,11);  u1,Q1 = SolveChannelModel(m1)
m2 = ChannelModel(1.,1.,3.,21);  u2,Q2 = SolveChannelModel(m2)
m3 = ChannelModel(1.,1.,3.,41);  u3,Q3 = SolveChannelModel(m3)
m4 = ChannelModel(1.,1.,3.,81);  u4,Q4 = SolveChannelModel(m4)
m5 = ChannelModel(1.,1.,3.,161); u5,Q5 = SolveChannelModel(m5)
m6 = ChannelModel(1.,1.,3.,321); u6,Q6 = SolveChannelModel(m6)

x2,y2 = ξη2xy(m2)
x4,y4 = ξη2xy(m4)
x6,y6 = ξη2xy(m6)

u6n  = reshape(u6,321,321)
u6_5 = u6n[1:2:end,1:2:end][:]
u6_4 = u6n[1:4:end,1:4:end][:]
u6_3 = u6n[1:8:end,1:8:end][:]
u6_2 = u6n[1:16:end,1:16:end][:]
u6_1 = u6n[1:32:end,1:32:end][:]

L2_1 = rmsd(u6_1,u1)
L2_2 = rmsd(u6_2,u2)
L2_3 = rmsd(u6_3,u3)
L2_4 = rmsd(u6_4,u4)
L2_5 = rmsd(u6_5,u5)

L∞_1 = maximum(abs.(u6_1.-u1))
L∞_2 = maximum(abs.(u6_2.-u2))
L∞_3 = maximum(abs.(u6_3.-u3))
L∞_4 = maximum(abs.(u6_4.-u4))
L∞_5 = maximum(abs.(u6_5.-u5))

Qd = abs.([Q1,Q2,Q3,Q4,Q5] .- Q6)
L2 = [L2_1,L2_2,L2_3,L2_4,L2_5]
L∞ = [L∞_1,L∞_2,L∞_3,L∞_4,L∞_5]

dlvl = [-10,-5,-2,-1,-0.5,-0.2,-0.1,0.1,0.2,0.5,1,2,5,10]

axs[1].tricontourf(x2[:],y2[:],u2,cmap="Blues",levels=(0:2:18)/100,extend="both")
axs[1].plot([m2.b/2,m2.b/2+m2.a],[0,m2.h],lw=0.5,c="k",linestyle="--")
axs[1].format(xlabel=L"x",ylabel=L"y",ltitle="(a) N = $(m2.n)",lrtitle="Q = $(@sprintf("%.5f",Q2))")

axs[2].tricontourf(x4[:],y4[:],u4,cmap="Blues",levels=(0:2:18)/100,extend="both")
axs[2].plot([m4.b/2,m4.b/2+m4.a],[0,m4.h],lw=0.5,c="k",linestyle="--")
axs[2].format(xlabel=L"x",ylabel=L"y",ltitle="(b) N = $(m4.n)",lrtitle="Q = $(@sprintf("%.5f",Q4))")

c1 = axs[3].tricontourf(x6[:],y6[:],u6,cmap="Blues",levels=(0:2:18)/100,extend="both")
axs[3].plot([m6.b/2,m6.b/2+m6.a],[0,m6.h],lw=0.5,c="k",linestyle="--")
axs[3].format(xlabel=L"x",ylabel=L"y",ltitle="(c) N = $(m6.n)",lrtitle="Q = $(@sprintf("%.5f",Q6))")

axs[4].tricontourf(x2[:],y2[:],(u2.-u6_2)*10000,cmap="RdBu",levels=dlvl,extend="both")
axs[4].plot([m2.b/2,m2.b/2+m2.a],[0,m2.h],lw=0.5,c="k",linestyle="--")
axs[4].format(xlabel=L"x",ylabel=L"y",ltitle="(d) N = $(m2.n)",lrtitle="|Q - Qr| = $(@sprintf("%.5f",Qd[2]))")

c2 = axs[5].tricontourf(x4[:],y4[:],(u4.-u6_4)*10000,cmap="RdBu",levels=dlvl,extend="both")
axs[5].plot([m4.b/2,m4.b/2+m4.a],[0,m4.h],lw=0.5,c="k",linestyle="--")
axs[5].format(xlabel=L"x",ylabel=L"y",ltitle="(e) N = $(m4.n)",lrtitle="|Q - Qr| = $(@sprintf("%.5f",Qd[4]))")

axs[6].plot(log10.([11,21,41,81,161]),log10.(Qd),label=L"|Q - \hat{Q}|",legend="ur")
axs[6].plot(log10.([11,21,41,81,161]),log10.(L2),label=L"L_2",legend="ur",legend_kw=Dict("frame"=>false,"ncol"=>1))
axs[6].plot(log10.([11,21,41,81,161]),log10.(L∞),label=L"L_\infty",legend="ur")
axs[6].plot([0,3],[-1.5,-7.5],lw=0.5,c="k",linestyle="--",label=L"|$\alpha$| = 2",legend="ur")
axs[6].format(xlim=(0,3),ylim=(-6,-3),xlabel="log N",ylabel=L"log $||u||_p$",ltitle="(f) Order of Convergence",suptitle="b = 1.0")

axs[3].colorbar(c1,loc="r",label="u")
axs[6].colorbar(c2,loc="r",label=L"u - $u_r$ / 10$^{-4}$")

fig.savefig(plotsdir("prj1_qn1f_b1d0.png"),transparent=false,dpi=250)