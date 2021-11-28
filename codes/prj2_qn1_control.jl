using DrWatson
@quickactivate "MIT6.339"
using MAT
using Printf
using PyCall, LaTeXStrings

pplt = pyimport("proplot")

include(srcdir("prj2","ThermalFin.jl"))

matdata = matread(srcdir("prj2_matlab","grids.mat"))
fine = matdata["fine"]
medi = matdata["medium"]
crse = matdata["coarse"]

pplt.close(); fig,axs = pplt.subplots([1,1,2],aspect=[8/4.4,1],axwidth=4,sharex=0,sharey=0)

u,Tfine = ThermalFin(fine,[0.4,0.6,0.8,1.2,0.1])
_,Tmed  = ThermalFin(medi,[0.4,0.6,0.8,1.2,0.1])
_,Tcrse = ThermalFin(crse,[0.4,0.6,0.8,1.2,0.1])

x = fine["coor"][:,1]
y = fine["coor"][:,2]
Ω1 = unique(Int.(fine["theta"][1])); nΩ1 = size(Ω1,1)
Ω2 = unique(Int.(fine["theta"][2])); nΩ2 = size(Ω2,1)
Ω3 = unique(Int.(fine["theta"][3])); nΩ3 = size(Ω3,1)
Ω4 = unique(Int.(fine["theta"][4])); nΩ4 = size(Ω4,1)
Ω0 = unique(Int.(fine["theta"][5])); nΩ0 = size(Ω0,1)
Γe = Int.(fine["theta"][6]); nΓe = size(Γe,1)
Γr = Int.(fine["theta"][7]); nΓr = size(Γr,1)
rot = [1,2,3,1]

c = axs[1].tricontourf(x[Ω1],y[Ω1],u[Ω1],levels=0:0.2:2,extend="both")
axs[1].tricontourf(x[Ω2],y[Ω2],u[Ω2],levels=0:0.2:2,extend="both")
axs[1].tricontourf(x[Ω3],y[Ω3],u[Ω3],levels=0:0.2:2,extend="both")
axs[1].tricontourf(x[Ω4],y[Ω4],u[Ω4],levels=0:0.2:2,extend="both")
axs[1].tricontourf(x[Ω0],y[Ω0],u[Ω0],levels=0:0.2:2,extend="both")
for ii = 1 : nΓe; axs[1].plot(x[Γe[ii,:]],y[Γe[ii,:]],c="k") end
                  axs[1].plot([-4,4],[0,0],c="brown")
for ii = 1 : nΓr; axs[1].plot(x[Γr[ii,:]],y[Γr[ii,:]],c="r") end
axs[1].format(xlim=(-4,4),ylim=(-0.2,4.2),ltitle="(a) Temperature Distribution (Medium Grid)")
axs[1].text(0.7,0.2,L"$T_{root}$ = " * @sprintf("%7.5f",Tmed))
axs[1].colorbar(c,loc="r")

axs[2].plot([-3,-1],[-7.5,-5.5],c="k",lw=0.5,linestyle="--")
axs[2].plot([-2.5,-1.5],[-7.5,-5.5],c="k",lw=0.5,linestyle="--")
axs[2].plot(log.([0.2,0.1]),log.(abs.([Tcrse,Tmed].-Tfine)))
axs[2].format(xlim=(-3,-1),ylim=(-7.5,-5.5),xlabel="log(h)",ylabel=L"log(e($T_{root}$)")
axs[2].format(ltitle=L"(b) Order of Convergence")

fig.savefig(plotsdir("prj2_qn1_control.png"),transparent=false,dpi=300)