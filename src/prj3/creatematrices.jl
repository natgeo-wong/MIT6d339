using DrWatson
@quickactivate "MIT6.339"
using LinearAlgebra

include(srcdir("prj3","quadrature.jl"))
include(srcdir("prj3","greensfunpot.jl"))

function ExteriorMatrices1(
    n :: Int
)

    A  = zeros(n,n)
    xc = zeros(n,2)
    xt = zeros(n,2)
    ψ  = zeros(n)
    θ0 = range(-π,π,length=n+1)
    θ1 = θ0[1:(end-1)]; θ2 = θ0[2:end]
    I  = -0.4237432928062354609011209
    σ  = zeros(n)

    for i = 1 : n
        A[i,i] -= π
        xc[i,1] = (cos(θ1[i]) + cos(θ2[i]))/2
        xc[i,2] = (sin(θ1[i]) + sin(θ2[i]))/2
        xt[i,1] = cos(θ1[i])
        xt[i,2] = sin(θ1[i])
        ψ[i]    = ∂u∂n1(atan(xc[i,2],xc[i,1]))
        σ[i]    = - (I/2 + ψ[i]) / π
    end

    for j = 1 : n, i = 1 : n
        A[i,j] -= 2*π/n * calcGreensFunction(xc,i,j)
    end

    return A,ψ,xc,xt,σ

end

function ExteriorMatrices2(
    n :: Int
)

    A  = zeros(n,n)
    xc = zeros(n,2)
    xt = zeros(n,2)
    ψ  = zeros(n)
    θ0 = range(-π,π,length=n+1)
    θ1 = θ0[1:(end-1)]; θ2 = θ0[2:end]

    for i = 1 : n
        A[i,i] -= π
        xc[i,1] = (cos(θ1[i]) + cos(θ2[i]))/2
        xc[i,2] = (sin(θ1[i]) + sin(θ2[i]))/2
        xt[i,1] = cos(θ1[i])
        xt[i,2] = sin(θ1[i])
        ψ[i]    = ∂u∂n2(xc[i,1],xc[i,2])
    end

    for j = 1 : n, i = 1 : n
        A[i,j] -= 2 * π / n * calcGreensFunction(xc,i,j)
    end

    return A,ψ,xc,xt

end

function ExteriorMatrices3(
    n :: Int;
    a :: Real
)

    A  = zeros(n,n)
    xc = zeros(n,2)
    xt = zeros(n,2)
    ψ  = zeros(n)
    θ0 = range(-π,π,length=n+1)
    θ1 = θ0[1:(end-1)]; θ2 = θ0[2:end]

    for i = 1 : n
        A[i,i] -= π
        xc[i,1] = (cos(θ1[i]) + cos(θ2[i])) / 2 * a
        xc[i,2] = (sin(θ1[i]) + sin(θ2[i])) / 2
        xt[i,1] = cos(θ1[i]) * a
        xt[i,2] = sin(θ1[i])
        ψ[i]    = ∂u∂n4(xc[i,1],xc[i,2],a)
    end

    for j = 1 : n, i = 1 : n
        A[i,j] -= 2 * π / n * calcGreensFunctionθ(xc,i,j,a) * 
                   sqrt((xc[i,1]/a)^2 + (a*xc[i,2])^2)
    end

    return A,ψ,xc,xt

end

function InteriorMatrices1(
    n :: Int
)

    A  = zeros(n,n)
    xc = zeros(n,2)
    xt = zeros(n,2)
    ψ  = zeros(n)
    θ0 = range(-π,π,length=n+1)
    θ1 = θ0[1:(end-1)]; θ2 = θ0[2:end]

    for i = 1 : n
        A[i,i] += π
        xc[i,1] = (cos(θ1[i]) + cos(θ2[i]))/2
        xc[i,2] = (sin(θ1[i]) + sin(θ2[i]))/2
        xt[i,1] = cos(θ1[i])
        xt[i,2] = sin(θ1[i])
        ψ[i]    = ∂u∂n3(xc[i,1],xc[i,2])
    end

    for j = 1 : n, i = 1 : n
        A[i,j] -= 2 * π / n * calcGreensFunction(xc,i,j)
    end

    for j = 1 : n
        A[n,j] = -2 * π / n * log(sqrt(xc[j,1]^2 + xc[j,2]^2))
    end

    ψ[n]   = 0

    return A,ψ,xc,xt

end