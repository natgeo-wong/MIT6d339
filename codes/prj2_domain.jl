using DrWatson
@quickactivate "MIT6.339"
using MAT
using PyCall, LaTeXStrings

pplt = pyimport("proplot")

matdata = matread(srcdir("prj2_matlab","grids.mat"))

pplt.close(); fig,axs = pplt.subplots(aspect=8/5,axwidth=3,ncols=2)

for imesh = 1 : 2
    
    if imesh == 1
        mesh = matdata["coarse"]
    elseif imesh == 2
        mesh = matdata["medium"]
    else
        mesh = matdata["fine"]
    end
    
    x = mesh["coor"][:,1]
    y = mesh["coor"][:,2]
    Ω1 = Int.(mesh["theta"][1]); nΩ1 = size(Ω1,1)
    Ω2 = Int.(mesh["theta"][2]); nΩ2 = size(Ω2,1)
    Ω3 = Int.(mesh["theta"][3]); nΩ3 = size(Ω3,1)
    Ω4 = Int.(mesh["theta"][4]); nΩ4 = size(Ω4,1)
    Ω0 = Int.(mesh["theta"][5]); nΩ0 = size(Ω0,1)
    Γe = Int.(mesh["theta"][6]); nΓe = size(Γe,1)
    Γr = Int.(mesh["theta"][7]); nΓr = size(Γr,1)
    rot = [1,2,3,1]

    for ii = 1 : nΩ1; axs[imesh].plot(x[Ω1[ii,rot]],y[Ω1[ii,rot]],lw=0.5,c="red3") end
    for ii = 1 : nΩ2; axs[imesh].plot(x[Ω2[ii,rot]],y[Ω2[ii,rot]],lw=0.5,c="orange3") end
    for ii = 1 : nΩ3; axs[imesh].plot(x[Ω3[ii,rot]],y[Ω3[ii,rot]],lw=0.5,c="yellow3") end
    for ii = 1 : nΩ4; axs[imesh].plot(x[Ω4[ii,rot]],y[Ω4[ii,rot]],lw=0.5,c="green3") end
    for ii = 1 : nΩ0; axs[imesh].plot(x[Ω0[ii,rot]],y[Ω0[ii,rot]],lw=0.5,c="gray7") end
    for ii = 1 : nΓe; axs[imesh].plot(x[Γe[ii,:]],y[Γe[ii,:]],c="k") end
                      axs[imesh].plot([-4,4],[0,0],c="brown")
    for ii = 1 : nΓr; axs[imesh].plot(x[Γr[ii,:]],y[Γr[ii,:]],c="r") end

    axs[imesh].format(xlim=(-4,4),ylim=(-0.2,4.8),suptitle="Finite Element Grid")

end

axs[1].format(ultitle="(a) Coarse",xlabel="x",ylabel="y")
axs[2].format(ultitle="(b) Medium")

fig.savefig(plotsdir("prj2_domain.png"),transparent=false,dpi=250)