using LinearAlgebra
using SparseArrays
using StatsBase

include(srcdir("prj1","qn2matrices.jl"))

abstract type SourceModel end

struct JacobiModel{FT} <: SourceModel
    N  :: Int
    S  :: Vector{Int}
    A  :: SparseMatrixCSC{Int,Int}
    R  :: SparseMatrixCSC{Float32,Int}
    f  :: SparseVector{Float32,Int}
    F  :: SparseVector{Float32,Int}
end

struct GaussSeidelModel{FT} <: SourceModel
    N  :: Int
    S  :: Vector{Int}
    A  :: SparseMatrixCSC{Int,Int}
    R  :: SparseMatrixCSC{Float32,Int}
    f  :: SparseVector{Float32,Int}
    F  :: Vector{Float32}
end

struct MultiGridModel{FT}
    N  :: Int
    S  :: Vector{Int}
    A  :: SparseMatrixCSC{Int,Int}
    f  :: SparseVector{Float32,Int}
    mR :: SourceModel
    mP :: SourceModel
    mC :: SourceModel
end

function JacobiModel(;
    N :: Int = 25,
    S :: Vector{Int} = [1,7,14,16],
    f  = 0,
    FT = Float64
)

    A  = createA(N)
    DI = createDinv(N)
    if iszero(f); f = createf(S,N) end
    R  = sparse(I,N^2,N^2) - DI * A
    F  = DI * f

    return JacobiModel{FT}(N,S,A,R,f,F)

end

function GaussSeidelModel(;
    N :: Int = 25,
    S :: Vector{Int} = [1,7,14,16],
    f  = 0,
    FT = Float64
)

    A  = createA(N)
    D  = Diagonal(A)
    U  = - UpperTriangular(A) + D
    L  = LowerTriangular(A)
    if iszero(f); f = createf(S,N) end
    R  = L \ U
    F  = L \ f

    return GaussSeidelModel{FT}(N,S,A,R,f,F)

end

function MultiGridModel(;
    N :: Int = 25,
    S :: Vector{Int} = [1,7,14,16],
    R :: String = "GS",
    P :: String = "GS",
    C :: String = "GS",
    f = 0,
    FT = Float64
)

    if R == "JC"
        mR = JacobiModel(N=N,S=S,f=f)
    else
        mR = GaussSeidelModel(N=N,S=S,f=f)
    end

    if P == "JC"
        mP = JacobiModel(N=N,S=S,f=f)
    else
        mP = GaussSeidelModel(N=N,S=S,f=f)
    end

    if C == "JC"
        mC = JacobiModel(N=N,S=S,f=f)
    else
        mC = GaussSeidelModel(N=N,S=S,f=f)
    end

    if iszero(f)
        f = createf(S,N)
    end

    return MultiGridModel{FT}(N,S,createA(N),f,mR,mP,mC)

end