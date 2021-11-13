using Dates
using LinearAlgebra
using Logging
using SparseArrays
using StatsBase

include(srcdir("prj1","qn2modelcreate.jl"))

function iterateModel(
    m :: SourceModel;
    niter :: Int,
    ω :: Real = 1
)

    ϕ = zeros(m.N^2)
    ϕr,_,_,_,_ = directsolveModel(m)
    ϵ = zeros(niter)
    iiter = 0
    ϵ[1] = rmsd(ϕ,ϕr)

    while iiter < niter
        ϕ = (m.R * ϕ + m.F) * ω + (1-ω) * ϕ
        iiter += 1
        ϵ[iiter] = rmsd(ϕ,ϕr)
    end

    ϕ = reshape(ϕ,m.N,m.N)
    lg = -1.5 * ϕ[1,:] + 2 * ϕ[2,:] - 0.5 * ϕ[3,:]
    rg = -0.5 * ϕ[end-2,:] + 2 * ϕ[end-1,:] -1.5 * ϕ[end,:]
    ug = -0.5 * ϕ[:,end-2] + 2 * ϕ[:,end-1] -1.5 * ϕ[:,end]
    bg = -1.5 * ϕ[:,1] + 2 * ϕ[:,2] - 0.5 * ϕ[:,3]

    return ϕ[:],ϵ,lg*(m.N-1),rg*(m.N-1),ug*(m.N-1),bg*(m.N-1)

end

function iterateModel(
    ϕ :: Array,
    m :: SourceModel;
    niter :: Int,
    ω :: Real = 1
)

    iiter = 0
    while iiter < niter
        ϕ = (m.R * ϕ + m.F) * ω + (1-ω) * ϕ
        iiter += 1
    end

    return ϕ

end

function directsolveModel(
    m :: SourceModel
)

    ϕ = m.A \ Vector(m.f)

    ϕ = reshape(ϕ,m.N,m.N)
    lg = -1.5 * ϕ[1,:] + 2 * ϕ[2,:] - 0.5 * ϕ[3,:]
    rg = -0.5 * ϕ[end-2,:] + 2 * ϕ[end-1,:] -1.5 * ϕ[end,:]
    ug = -0.5 * ϕ[:,end-2] + 2 * ϕ[:,end-1] -1.5 * ϕ[:,end]
    bg = -1.5 * ϕ[:,1] + 2 * ϕ[:,2] - 0.5 * ϕ[:,3]

    return ϕ[:],lg*(m.N-1),rg*(m.N-1),ug*(m.N-1),bg*(m.N-1)

end

function iterateMultiGrid(
    m  :: MultiGridModel;
    ω  :: Real,
    ν1 :: Int = 1,
    ν2 :: Int = 1,
    νc :: Int = ν1+ν2,
    nν :: Int,
    niter :: Int,
)

    ϕ = zeros(m.N^2); ϕr,_,_,_,_ = directsolveModel(m.mR)
    ϵ = zeros(niter)

    iiter = 0
    while iiter < niter
        iν = 1
        ϕ = recursiveMultiGrid(Vector(ϕ),m,ω,ν1,ν2,νc,iν,nν)
        iiter += 1; ϵ[iiter] = rmsd(ϕ,ϕr)
        if iszero(mod(iiter,50))
            @info "$(now()) - MultiGrid iteration number $iiter"
        end
    end

    return ϕ,ϵ
    
end

function recursiveMultiGrid(
    ϕ  :: Array,
    m  :: MultiGridModel,
    ω  :: Real,
    ν1 :: Int,
    ν2 :: Int,
    νc :: Int,
    iν :: Int,
    nν :: Int,
)

    S = m.S
    ϕ = iterateModel(ϕ,m.mR,niter=ν1,ω=ω)
    rA,rN = restrictionMatrix(m.N)
    resh  = m.f - m.A * ϕ
    res2h = rA * resh
    newm  = MultiGridModel(N=rN,S=S,f=res2h)
    
    if (iν < nν-1) && !iszero(mod((newm.N+1)/2,2))
        iν += 1
        e2h = recursiveMultiGrid(zeros(rN^2),newm,ω,ν1,ν2,νc,iν,nν)
    else
        if !iszero(νc)
            e2h = iterateModel(zeros(rN^2),newm.mC,niter=νc,ω=ω)
        else
            e2h = newm.A \ Vector(newm.f);
        end
    end

    pA,_ = prolongateMatrix(newm.N)
    eh  = pA * e2h
    ϕ   = ϕ + Vector(eh)
    ϕ = iterateModel(Vector(ϕ),m.mP,niter=ν2,ω=ω)

    return ϕ

end

function ϕ2∂ₙϕ(
    ϕ :: Array
)

    N = Int(sqrt(length(ϕ)))
    ϕ = reshape(ϕ,N,N)
    lg = -1.5 * ϕ[1,:] + 2 * ϕ[2,:] - 0.5 * ϕ[3,:]
    rg = -0.5 * ϕ[end-2,:] + 2 * ϕ[end-1,:] -1.5 * ϕ[end,:]
    ug = -0.5 * ϕ[:,end-2] + 2 * ϕ[:,end-1] -1.5 * ϕ[:,end]
    bg = -1.5 * ϕ[:,1] + 2 * ϕ[:,2] - 0.5 * ϕ[:,3]

    return lg*(N-1),rg*(N-1),ug*(N-1),bg*(N-1)

end