using DrWatson
@quickactivate "MIT6.339"
using PrettyTables
using PyCall, LaTeXStrings
pplt = pyimport("proplot")

include(srcdir("prj1qn2functions.jl"))

m = GaussSeidelModel(N=25)
ϕ,ϵ,lg,rg,ug,bg = iterateModel(m,niter=750,ω=1.7)

fmat = cat(collect(1:25),lg,rg,bg,ug,dims=2)

head = ["Node","Left","Right","Bottom","Top"];

pretty_table(
    fmat,head,
    alignment=[:c,:c,:c,:c,:c],
    crop = :none, tf = tf_compact
);