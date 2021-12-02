using DrWatson
@quickactivate "MIT6.339"
using LinearAlgebra

include(srcdir("prj3","quadrature.jl"))

function calcGreensFunction(
    xvec :: Array{<:Real,2},
    i :: Int,  j :: Int
)

    if i != j

        xi = xvec[i,1]; yi = xvec[i,2]
        xj = xvec[j,1]; yj = xvec[j,2]

        return (xi * (xi-xj) + yi * (yi-yj)) / ((xi-xj)^2 + (yi-yj)^2)

    else

        return 0.5

    end

end

function calcGreensFunctionθ(
    xvec :: Array{<:Real,2},
    i :: Int,  j :: Int, a :: Real
)

    xi = xvec[i,1]; yi = xvec[i,2]

    if i != j

        xj = xvec[j,1]; yj = xvec[j,2]
        nx = xi / sqrt(xi^2/a^2 + a^2*yi^2) / a
        ny = yi / sqrt(xi^2/a^2 + a^2*yi^2) * a

        return (nx * (xi-xj) + ny * (yi-yj)) / ((xi-xj)^2 + (yi-yj)^2)

    else

        return 0.5 * a / sqrt(a^2*yi^2+xi^2/a^2)^3

    end

end

function ∂u∂n1(θ)

    return 1 / (3+2*cos(θ)+cos(2*θ))

end

function ∂u∂n2(x,y)

    x2  = x^2
    xy  = sqrt(x2+y^2)
    xym = x2+(y-0.5)^2
    xyp = x2+(y+0.5)^2
    ∂u∂x = x2  / xy * (1/(xyp) - 1/(xym))
    ∂u∂y = y/2 / xy * ((2*y+1)/(xyp) - (2*y-1)/(xym))

    return ∂u∂x + ∂u∂y

end

function ∂u∂n3(x,y)

    x2  = x^2
    xy  = sqrt(x2+y^2)
    xym = x2+(y-2)^2
    xyp = x2+(y+2)^2
    ∂u∂x = x2 / xy * (1/(xyp) - 1/(xym))
    ∂u∂y = y  / xy * ((y+2)/(xyp) - (y-2)/(xym))

    return ∂u∂x + ∂u∂y

end

function ∂u∂n4(x,y,a)

    x2  = x^2/a^2
    xy  = sqrt(x2+y^2*a^2)
    xym = x2*a^2+(y-0.5)^2
    xyp = x2*a^2+(y+0.5)^2
    ∂u∂x = x2*a  / xy * (1/(xyp) - 1/(xym))
    ∂u∂y = y*a/2 / xy * ((2*y+1)/(xyp) - (2*y-1)/(xym))

    return ∂u∂x + ∂u∂y

end

function upotnumeric(
    σ  :: Vector{<:Real},
    xc :: Array{<:Real,2},
    xj :: Vector{<:Real}
)

    n  = length(σ)
    xd = xc[:,1] .- xj[1]
    yd = xc[:,2] .- xj[2]
    dd = log.(sqrt.(xd.^2 .+ yd.^2))
    up = - 2 * π / n * sum(σ .* dd)

    return up

end

function upotnumericθ(
    σ  :: Vector{<:Real},
    xc :: Array{<:Real,2},
    xj :: Vector{<:Real};
    a  :: Real
)

    n  = length(σ)
    xd = xc[:,1] .- xj[1]
    yd = xc[:,2] .- xj[2]
    dd = log.(sqrt.(xd.^2 .+ yd.^2))
    up = - 2 * π / n * sum(σ .* dd .* sqrt.(xc[:,2].^2*a^2 .+ xc[:,1].^2/a^2))

    return up

end

function upotanalyticqn1(xj :: Vector{<:Real})

    x = xj[1]
    y = xj[2]
    return (log(x^2+(y+0.5)^2) - log(x^2+(y-0.5)^2)) * 0.5

end

function upotanalyticqn2(xj :: Vector{<:Real})

    x = xj[1]
    y = xj[2]
    return (log(x^2+(y+2)^2) - log(x^2+(y-2)^2)) * 0.5

end