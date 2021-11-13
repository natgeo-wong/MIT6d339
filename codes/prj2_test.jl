using DrWatson
@quickactivate "MIT6.339"
using MAT
using PyCall, LaTeXStrings

pplt = pyimport("proplot")

include(srcdir("prj2","ThermalFin.jl"))

matdata = matread(srcdir("prj2_matlab","grids.mat"))
mesh = matdata["medium"]

pplt.close(); fig,axs = pplt.subplots(aspect=8/5,axwidth=3)

u = ThermalFin(mesh,[0.4,0.6,0.8,1.2,0.1])

x = mesh["coor"][:,1]
y = mesh["coor"][:,2]
Ω1 = unique(Int.(mesh["theta"][1])); nΩ1 = size(Ω1,1)
Ω2 = unique(Int.(mesh["theta"][2])); nΩ2 = size(Ω2,1)
Ω3 = unique(Int.(mesh["theta"][3])); nΩ3 = size(Ω3,1)
Ω4 = unique(Int.(mesh["theta"][4])); nΩ4 = size(Ω4,1)
Ω0 = unique(Int.(mesh["theta"][5])); nΩ0 = size(Ω0,1)
Γe = Int.(mesh["theta"][6]); nΓe = size(Γe,1)
Γr = Int.(mesh["theta"][7]); nΓr = size(Γr,1)
rot = [1,2,3,1]

c = axs[1].tricontourf(x[Ω1],y[Ω1],u[Ω1],levels=0:0.1:1.5,extend="both")
axs[1].tricontourf(x[Ω2],y[Ω2],u[Ω2],levels=0:0.1:1.5,extend="both")
axs[1].tricontourf(x[Ω3],y[Ω3],u[Ω3],levels=0:0.1:1.5,extend="both")
axs[1].tricontourf(x[Ω4],y[Ω4],u[Ω4],levels=0:0.1:1.5,extend="both")
axs[1].tricontourf(x[Ω0],y[Ω0],u[Ω0],levels=0:0.1:1.5,extend="both")
for ii = 1 : nΓe; axs[1].plot(x[Γe[ii,:]],y[Γe[ii,:]],c="k") end
                  axs[1].plot([-4,4],[0,0],c="brown")
for ii = 1 : nΓr; axs[1].plot(x[Γr[ii,:]],y[Γr[ii,:]],c="r") end
axs[1].format(xlim=(-4,4),ylim=(-0.5,4.5))

fig.colorbar(c,loc="r")
fig.savefig(plotsdir("prj2_test.png"),transparent=false,dpi=150)