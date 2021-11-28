using LinearAlgebra
using SparseArrays

function ThermalFin(mesh, mu)
    #
    # -------------------------------------------------------------------------
    #   Computes the temperature distribution and root temperature for a fin
    #   using the Finite Element Method.
    # -------------------------------------------------------------------------
    #
    #   INPUT   mesh    grid lebel (coarse, medium, or fine)
    #           mu      thermal conductivities of sections and Biot number 1x5
    #
    #   OUTPUT  u       temperature disctribution in the fin
    #           Troot   root temperature
    #
    # -------------------------------------------------------------------------
    
    
    # Parameters setup: rearranges the values in mu
    kappa = ones(6)
    kappa[1:4] .= mu[1:4]
    kappa[6] = mu[5]

    nnodes = Int.(mesh["nodes"])

    # Initialization
    A = spzeros(nnodes,nnodes)
    F = zeros(nnodes)
    L = zeros(1,nnodes)
    A3by3 = zeros(3,3)
    A2by2 = ones(2,2) + I(2)
    x3by3 = ones(3,3)
    
    
    # Domain Interior
    for i = 1 : 5, n = 1 : size(mesh["theta"][i],1)
        phi = Int.(mesh["theta"][i][n,:])
        
        for j = 1 : 3
            x3by3[j,2:3] .= mesh["coor"][phi[j],:]
        end

        x1k = x3by3[1,2]; x2k = x3by3[2,2]; x3k = x3by3[3,2]
        y1k = x3by3[1,3]; y2k = x3by3[2,3]; y3k = x3by3[3,3]

        ThkA = 0.5 * abs(x1k*(y2k-y3k) + x2k*(y3k-y1k) + x3k*(y1k-y2k))

        cmat = inv(x3by3)

        for β = 1 : 3, α = 1 : 3
            A3by3[α,β] = ThkA * kappa[i] * (cmat[2,α] * cmat[2,β] + cmat[3,α] * cmat[3,β])
        end

        for β = 1 : 3, α = 1 : 3
            A[phi[α],phi[β]] = A[phi[α],phi[β]] + A3by3[α,β]
        end
    end
    
    
    # Boundaries (not root)
    for n = 1 : size(mesh["theta"][6],1)
        phi = Int.(mesh["theta"][6][n,:])
        
        x1k = mesh["coor"][phi[1],1]; x2k = mesh["coor"][phi[2],1]
        y1k = mesh["coor"][phi[1],2]; y2k = mesh["coor"][phi[2],2]
        xyL2 = sqrt((y2k-y1k)^2 + (x2k-x1k)^2)
    
        for β = 1 : 2, α = 1 : 2
            A[phi[α],phi[β]] = A[phi[α],phi[β]] + A2by2[α,β] * xyL2 * kappa[6] / 6
        end
    end
    
    
    # Root Boundary
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
    
    
    # Solve for the temperature distribution
    u = A \ F
    Troot = L * u

    return u, Troot[1]
    
end