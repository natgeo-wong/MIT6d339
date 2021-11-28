using DrWatson
@quickactivate "MIT6.339"
using JLD2
using Printf
using PyCall, LaTeXStrings

include(srcdir("prj2","ReducedBaseOnline.jl"))

@load "ReducedBasis.jld2" N ANq FN LN

_,TN_1 = ReducedBasisOnline([0.4,0.6,0.8,1.2,0.1],N,ANq,FN,LN)
_,TN_2 = ReducedBasisOnline([1.8,4.2,5.7,2.9,0.3],N,ANq,FN,LN)