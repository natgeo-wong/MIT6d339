function ReducedBasisOnline(mu,N,ANq,FN,LN)
#
# -------------------------------------------------------------------------
#
#   ReducedBasisOnline.m
#
# -------------------------------------------------------------------------
#   Reduced basis on-line: computes the temperature distribution and root
#   temperature using the reduced basis computed by the off-line part.
#
#   The reduced basis data ReducedBasis.mat should have been previsously
#   loaded using the comand "load ReducedBasis".
# -------------------------------------------------------------------------
#
#   INPUT   mu      thermal conductivities of sections and Biot number 1x5
#           N       dimension of the reduced basis
#           ANq     reduced matrix
#           FN      reduced source
#
#   OUTPUT  uN      temperature disctribution in the fin from reduced basis
#           TrootN  root temperature from reduced basis
#
# -------------------------------------------------------------------------


    # Parameters - (domain 5 is the root)
    sigma = ones(6);
    sigma[1:4] = mu[1:4];    # Bi = mu(6)
    sigma[6] = mu[5];

    # Construct of A_N
    AN = zeros(N,N)
    for iQ = 1 : length(ANq)
        AN += ANq[iQ] * sigma[iQ]
    end

    # Solve the system
    uN = AN \ FN;

    # Compute the temperature at the root
    TrootN = transpose(LN) * uN

    return uN, TrootN

end