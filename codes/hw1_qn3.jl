using DrWatson
@quickactivate "MIT6.339"
using Distributions

using PyCall, LaTeXStrings
pplt = pyimport("proplot")

curveinitial(x) = 20 / sqrt(2*π) * exp(-200 * (x-0.2)^2)
curvefinal(x,N,M,p) = 1 / sqrt(4*π*(N/2/M^2+0.00125)) * exp(-(x - 0.2 - (2p-1)*N/M)^2/4/(N/2/M^2+0.00125))

## Random Walk Model for Convection-Diffusion on a [0,1] periodic domain

## Parameters
M = 400     # Number of subdivisions in x
N = 300     # Number of timesteps
p = 0.8     # Probability of moving to the right

X = 0 : (1/M) : (1 - 1/M)

d = Normal(0.2,0.05)
dcdf = cdf.(d,0 : (1/M) : 1)

PN = M*(dcdf[2:end].-dcdf[1:(end-1)]) # centered at 0.2 and with std 0.05
PN1 = zeros(length(PN))

pplt.close(); f,a = pplt.subplots(aspect=2,axwidth=4)

a[1].plot(X,PN,label="Initial State",legend="ur",legend_kw=Dict("frame"=>false,"ncol"=>1))
a[1].plot(X,curveinitial.(X),c="k",linestyle="--")

for ii = 1 : N

    for iM = 2 : (M-1)
        PN1[iM] = p * PN[iM-1] + (1-p) * PN[iM+1]
    end
    PN1[1] =  p*PN[end] + (1-p)*PN[2]
    PN1[end] =  p*PN[end-1] + (1-p)*PN[1]

    for iM = 1 : M
        PN[iM] = PN1[iM]
    end

end

a[1].plot(X,PN,legend="ur",label="Final State")
a[1].plot(X,curvefinal.(X,N,M,p),c="k",linestyle="--",legend="ur",label="Analytical Solutions")
a[1].format(xlim=(0,1),ylim=(0,10),xlabel="X",ylabel="Y")

f.savefig(plotsdir("hw1_qn3.png"),transparent=false,dpi=150)