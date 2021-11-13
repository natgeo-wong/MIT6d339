using LinearAlgebra
using NumericalIntegration
using SparseArrays

struct ChannelModel{FT<:Real}
    n  :: Int
    a  :: FT
    b  :: FT
    d  :: FT
    h  :: FT
    l  :: FT
    c1 :: Vector{FT}
    c2 :: Array{FT,2}
    c3 :: Vector{FT}
    c5 :: Vector{FT}
    J  :: Vector{FT}
    ξ  :: Vector{FT}
    η  :: Vector{FT}
end

function ChannelModel(b::Real,h::Real,l::Real,n::Int;FT = Float64)

    d = 0.5 * (l-b)
    a = sqrt(0.25*(l-b)^2 - h^2)
    ξ = collect(0:n-1) / (n-1)
    η = collect(0:n-1) / (n-1)

    J  = b * h / 2 .+ η * h * a
    c1 = ξ.^2 * a^2 .+ h^2
    c2 = b * ξ * a / 2 .+ ξ .* η' * a^2
    c3 = (b/2 .+ a*η).^2
    c5 = 2 * a^2 * ξ

    return ChannelModel{FT}(n,a,b,d,h,l,c1,c2,c3,c5,J,ξ,η)

end

function SolveChannelModel(m::ChannelModel)

    A = createA(m)
    f = createf(m)
    u = A \ f
    Q = zeros(m.n^2)
    Q2 = zeros(m.n-1,m.n-1)
    n2 = (m.n-1)^2

    for j = 1 : m.n, i = 1 : m.n
        idiag = i + (j-1) * m.n
        Q[idiag] = u[idiag] * m.h * (m.b/2 + m.a * m.η[j]) / n2
    end

    Q = reshape(Q,m.n,m.n)
    for j = 1 : (m.n-1), i = 1 : (m.n-1)

        Q2[i,j] = (Q[i,j] + Q[i,j+1] + Q[i+1,j] + Q[i+1,j+1])/4

    end

    return u,2*sum(Q2)

end

function ξη2xy(m::ChannelModel)

    ξ = m.ξ
    η = m.η

    x = ξ .* (m.b/2 .+ m.a*η')
    y = zeros(size(x))
    y .= η' * m.h

    return x,y

end

function createA(m::ChannelModel)

    A = spzeros(m.n^2,m.n^2) .+ sparse(I,m.n^2,m.n^2)
    n1 = m.n-1
    n2 = (m.n-1)^2
    xm = m.a + m.b / 2

    for j = 2 : (m.n-1), i = 2 : (m.n-1)

        idiag = i + (j-1) * m.n
        A[idiag,idiag]   = -2 * (m.c1[i] + m.c3[j]) * n2
        A[idiag,idiag-1] = m.c1[i] * n2 - m.c5[i] * n1 / 2
        A[idiag,idiag+1] = m.c1[i] * n2 + m.c5[i] * n1 / 2
        A[idiag,idiag-m.n] = m.c3[j] * n2
        A[idiag,idiag+m.n] = m.c3[j] * n2
        A[idiag,idiag-m.n-1] = -0.5 * m.c2[i,j] * n2
        A[idiag,idiag+m.n+1] = -0.5 * m.c2[i,j] * n2
        A[idiag,idiag+m.n-1] =  0.5 * m.c2[i,j] * n2
        A[idiag,idiag-m.n+1] =  0.5 * m.c2[i,j] * n2

    end

    for j = 2 : (m.n-1)

        idiag = 1 + (j-1) * m.n
        A[idiag,idiag]   = 3
        A[idiag,idiag+1] = -4
        A[idiag,idiag+2] = 1

    end

    for i = 1 : (m.n-2)

        idiag = i + (m.n-1) * m.n
        A[idiag,idiag]       =  3 * xm + (i-1) * m.a / n1 * 3
        A[idiag,idiag-m.n]   = -4 * xm
        A[idiag,idiag-2*m.n] =  1 * xm
        A[idiag,idiag+1]     = (i-1) * m.a / n1 * -4
        A[idiag,idiag+2]     = (i-1) * m.a / n1 * 1

    end

    for i = (m.n-1) : (m.n-1)

        idiag = i + (m.n-1) * m.n
        A[idiag,idiag]       =  3 * xm
        A[idiag,idiag-m.n]   = -4 * xm
        A[idiag,idiag-2*m.n] =  1 * xm
        A[idiag,idiag-1]     = (i-1) * m.a / n1

    end

    return A

end

function createf(m)

    f = zeros(m.n^2)

    for j = 2 : (m.n-1), i = 2 : (m.n-1)

        idiag = i + (j-1) * m.n
        f[idiag] = - m.J[j]^2

    end

    return f

end