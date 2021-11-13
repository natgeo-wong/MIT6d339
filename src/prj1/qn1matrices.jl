using SparseArrays

include(srcdir("prj1","qn1models.jl"))

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