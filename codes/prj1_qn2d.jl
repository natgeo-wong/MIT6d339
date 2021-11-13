using DrWatson
@quickactivate "MIT6.339"
using MAT
using PrettyTables
using StatsBase

using PyCall, LaTeXStrings
pplt = pyimport("proplot")

include(srcdir("prj1qn2functions.jl"))

fluxr = matread(srcdir("flux_unknown.mat"))["flux"]
i = 0
ϵ = zeros(1820)
smat = zeros(4,1820)

for is1 = 1 : 13, is2 = (is1+1) : 14, is3 = (is2+1) : 15, is4 = (is3+1) : 16

    global i += 1
    m = MultiGridModel(N=25,S=[is1,is2,is3,is4])
    ϕ,_ = iterateMultiGrid(m,niter=20,ω=0.8,ν1=1,ν2=1,νc=50,nν=2)
    lgn,rgn,ugn,bgn = ϕ2∂ₙϕ(ϕ)
    fluxn = cat(lgn,rgn,bgn,ugn,dims=2)
    ϵ[i] = rmsd(fluxn,fluxr)
    smat[:,i] = [is1,is2,is3,is4]
    @info i

end

mf = GaussSeidelModel(N=25,S=Int.(smat[:,argmin(ϵ)]))
ϕf,lgf,rgf,ugf,bgf = directsolveModel(mf)
ϕf = reshape(ϕf,mf.N,mf.N)

x = (0 : (mf.N-1)) / (mf.N-1)

pplt.close(); fig,axs = pplt.subplots(axwidth=2)

c = axs[1].contourf(x,x,ϕf'*1000,levels=(2:10)*2.5,extend="both")
axs[1].plot([1,2,2,1,1]/6,[1,1,2,2,1]/6,c="k")
axs[1].plot([2,3,3,2,2]/6,[2,2,3,3,2]/6,c="k")
axs[1].plot([3,4,4,3,3]/6,[3,3,4,4,3]/6,c="k")
axs[1].plot([4,5,5,4,4]/6,[2,2,3,3,2]/6,c="k")
axs[1].format(xlocator=0:1/6:1,xticklabels=[],ylocator=0:1/6:1,yticklabels=[])
axs[1].format(suptitle="Unknown Flux Solution")
axs[1].colorbar(c,loc="r",locator=5:5:25,label=L"$\phi$ / 10$^{-3}$")

fig.savefig(plotsdir("prj1_qn2d.png"),transparent=false,dpi=250)

head = ["Node","Left","Right","Bottom","Top"];

pretty_table(
    cat(collect(1:25),lgf,rgf,bgf,ugf,dims=2),head,
    alignment=[:c,:c,:c,:c,:c],
    crop = :none, tf = tf_compact
);