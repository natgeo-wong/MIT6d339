using SparseArrays

function createA(N)

    N2 = (N-1)^2
    A  = spzeros(N^2,N^2) .+ Diagonal(ones(N^2)) * N2

    for j = 2 : (N-1), i = 2 : (N-1)

        idiag = i + (j-1) * N
        A[idiag,idiag]   =  N2 * 4
        A[idiag,idiag-1] = -N2
        A[idiag,idiag+1] = -N2
        A[idiag,idiag-N] = -N2
        A[idiag,idiag+N] = -N2

    end

    return A

end

function createf(
    S :: Vector{Int},
    N :: Int,
)

    n = Int((N-1) / 6)
    f = spzeros(N^2)

    for iS in S

        s1 = mod(iS,4); if iszero(s1); s1 = 4 end
        s2 = Int((iS - s1)/4 + 1)

        for j = (n*s1+1) : (n*(s1+1)+1), i = (n*s2+1) : (n*(s2+1)+1)

            f[i+(j-1)*N] = 1

        end

    end

    return f

end

function createDinv(N)

    N2   = (N-1)^2
    Dinv = spzeros(N^2,N^2) .+ Diagonal(ones(N^2)) / N2

    for j = 2 : (N-1), i = 2 : (N-1)

        idiag = i + (j-1) * N
        Dinv[idiag,idiag] = 0.25 / N2

    end

    return Dinv

end

function restrictionMatrix(N :: Int, factor :: Int = 2)

    if !iszero(mod((N-1)/factor,factor))
        error("N cannot be subdivided from a h-grid to a $(factor)h-grid.")
    end
    
    nN = Int((N-1)/factor + 1)
    A  = spzeros(nN*nN,N*N)

    for j = 1 : nN, i = 1 : nN

        ndiag = i + (j-1) * nN
        odiag = (i-1) * factor + 1 + (j-1) * factor * N

        A[ndiag,odiag] = 1

    end

    return A,nN

end

function prolongateMatrix(N :: Int)

    nN = (N-1)*2 + 1
    A  = spzeros(nN*nN,N*N)

    for j = 1 : N, i = 1 : N

        idiag1 = (2*i - 1) + (2*j - 2) * nN
        A[idiag1,i+(j-1)*N] = 1

        if i != N
            idiag2 = 2*i + (2*j - 2) * nN
            A[idiag2,i+(j-1)*N] = 0.5
            A[idiag2,i+1+(j-1)*N] = 0.5
        end

        if j != N
            idiag3 = (2*i - 1) + (2*j - 1) * nN
            A[idiag3,i+(j-1)*N] = 0.5
            A[idiag3,i+j*N] = 0.5
        end

        if (i != N) && (j != N)
            idiag4 = 2*i + (2*j - 1) * nN
            A[idiag4,i+(j-1)*N] = 0.25
            A[idiag4,i+j*N] = 0.25
            A[idiag4,i+1+(j-1)*N] = 0.25
            A[idiag4,i+1+j*N] = 0.25
        end

    end

    return A,nN

end