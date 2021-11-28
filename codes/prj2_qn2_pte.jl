using DrWatson
@quickactivate "MIT6.339"
using JLD2
using Optim
using Printf
using PyCall, LaTeXStrings

pplt = pyimport("proplot")

include(srcdir("prj2","ReducedBaseOnline.jl"))

@load "ReducedBasis.jld2" N ANq FN LN

function Tcost(Bi)
    _,T = ReducedBasisOnline([0.4,0.6,0.8,1.2,Bi[1]],N,ANq,FN,LN)

    return Bi[1] * 0.2 + T
end
fout = optimize(Tcost,[2.0],Newton())
fout.minimizer[1]

Bivec = 0.01:0.01:10
Troot = zeros(length(Bivec))
for iBi = 1 : length(Bivec)
    _,Troot[iBi] = ReducedBasisOnline([0.4,0.6,0.8,1.2,Bivec[iBi]],N,ANq,FN,LN)
end
cost  = 0.2 * Bivec .+ Troot

_,Tcostmin = ReducedBasisOnline([0.4,0.6,0.8,1.2,fout.minimizer[1]],N,ANq,FN,LN)
costmin = fout.minimizer[1] * 0.2 + Tcostmin

pplt.close(); fig,axs = pplt.subplots()

axs[1].plot(Bivec,cost,label="Cost",legend="r",legend_kw=Dict("frame"=>false,"ncol"=>1))
axs[1].plot(Bivec,Troot,label="Biot Number",legend="r")
axs[1].scatter(fout.minimizer[1],costmin)
axs[1].scatter(fout.minimizer[1],Tcostmin)

axs[1].plot([fout.minimizer[1],fout.minimizer[1]],[0,4],c="k",linestyle="--",lw=0.5)
axs[1].plot([0,10],[Tcostmin,Tcostmin],c="k",linestyle="--",lw=0.5)
axs[1].plot([0,10],[costmin,costmin],c="k",linestyle="--",lw=0.5)

axs[1].text(fout.minimizer[1]+0.3,3.7,"$(@sprintf("Bi = %7.5f",fout.minimizer[1]))")
axs[1].text(7,costmin+0.07,"$(@sprintf("C = %7.5f",costmin))")
axs[1].text(7,Tcostmin-0.22,"$(@sprintf("T = %7.5f",Tcostmin))")

axs[1].format(xlim=(0,10),ylim=(0,4),xlabel="Bi")

fig.savefig(plotsdir("prj2_qn2_pte.png"),transparent=false,dpi=150)