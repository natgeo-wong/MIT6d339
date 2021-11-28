using DrWatson
@quickactivate "MIT6.339"
using DelimitedFiles
using JLD2
using MAT
using SparseArrays

# -------------------------------------------------------------------------
#
#   ReducedBasisOffline.m
#
# -------------------------------------------------------------------------
#   Reduced basis off-line: computes the reduced basis matrices and source 
#   vector from a sample.
# -------------------------------------------------------------------------

include(srcdir("prj2","ThermalFin.jl"))

# Load grid triangulations
mesh = matread(srcdir("prj2_matlab","grids.mat"))["coarse"]
μmat = readdlm(srcdir("prj2_matlab","sn.dat"))
nμ   = size(μmat,1)
nnodes = Int.(mesh["nodes"])

# Solve the sample cases and construct Z
Z = zeros(nnodes,N);
F = zeros(nnodes)
L = zeros(nnodes)

for iμ = 1 : nμ
    Z[:,iμ],_ = ThermalFin(mesh,μmat[iμ,:])
end

for n = 1 : size(mesh["theta"][7],1)
    phi = Int.(mesh["theta"][7][n,:])
    
    x1k = mesh["coor"][phi[1],1]; x2k = mesh["coor"][phi[2],1]
    y1k = mesh["coor"][phi[1],2]; y2k = mesh["coor"][phi[2],2]
    xyL2 = sqrt((y2k-y1k)^2 + (x2k-x1k)^2)

    for ii = 1 : 2
        F[phi[ii]] = F[phi[ii]] + xyL2 / 2
        L[phi[ii]] = L[phi[ii]] + 0.5 * xyL2
    end
end

# Initialization of reduced-base matrix and vector
ANq = [zeros(N,N),zeros(N,N),zeros(N,N),zeros(N,N),zeros(N,N),zeros(N,N)]
A3by3 = zeros(3,3)
A2by2 = ones(2,2) + I(2)
x3by3 = ones(3,3)

for iQ = 1 : 5

    A = spzeros(nnodes,nnodes)
    
    for n = 1 : size(mesh["theta"][iQ],1)

        phi = Int.(mesh["theta"][iQ][n,:])
        
        for j = 1 : 3
            x3by3[j,2:3] .= mesh["coor"][phi[j],:]
        end

        x1k = x3by3[1,2]; x2k = x3by3[2,2]; x3k = x3by3[3,2]
        y1k = x3by3[1,3]; y2k = x3by3[2,3]; y3k = x3by3[3,3]

        ThkA = 0.5 * abs(x1k*(y2k-y3k) + x2k*(y3k-y1k) + x3k*(y1k-y2k))

        cmat = inv(x3by3)

        for β = 1 : 3, α = 1 : 3
            A3by3[α,β] = ThkA * (cmat[2,α] * cmat[2,β] + cmat[3,α] * cmat[3,β])
        end

        for β = 1 : 3, α = 1 : 3
            A[phi[α],phi[β]] = A[phi[α],phi[β]] + A3by3[α,β]
        end

    end

    ANq[iQ] = transpose(Z) * A * Z

end

A6 = spzeros(nnodes,nnodes)

for n = 1 : size(mesh["theta"][6],1)
    phi = Int.(mesh["theta"][6][n,:])
    
    x1k = mesh["coor"][phi[1],1]; x2k = mesh["coor"][phi[2],1]
    y1k = mesh["coor"][phi[1],2]; y2k = mesh["coor"][phi[2],2]
    xyL2 = sqrt((y2k-y1k)^2 + (x2k-x1k)^2)

    for β = 1 : 2, α = 1 : 2
        A6[phi[α],phi[β]] = A6[phi[α],phi[β]] + A2by2[α,β] * xyL2 / 6
    end
end

ANq[6] = transpose(Z) * A6 * Z

# Define FN
FN = zeros(N); FN = transpose(Z) * F

# Define LN
LN = zeros(N); LN = transpose(Z) * L

# Save the reduced matrix and vector
@save "ReducedBasis.jld2" N ANq FN LN