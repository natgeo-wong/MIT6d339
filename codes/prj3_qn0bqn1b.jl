using DrWatson
@quickactivate "MIT6.339"
using SpecialFunctions

using PyCall, LaTeXStrings
pplt = pyimport("proplot")

include(srcdir("prj3","quadrature.jl"))

f1(x) = 3*x
f2(x) = sin(x)
f3(x) = exp(2*cos(2*π*x))
f4(x) = abs(cos(x))

g1(x) = 3*x^2/2; g1a = g1(10) - g1(2)
g2(x) = -cos(x); g2a = g2(π) - g2(0)
                 g3a = besseli(0,2)
g4(x) = sin(x);  g4a = 4 * (g4(π/2) - g4(0))

n = 1 : 10000

f1a = zeros(length(n))
f2a = zeros(length(n))
f3a = zeros(length(n))
f4a = zeros(length(n))

for ii = 1 : length(n)

    f1a[ii] = abs(comptrap(f1,a=2,b=10,n=n[ii]) - g1a)
    f2a[ii] = abs(comptrap(f2,a=0,b=π,n=n[ii]) - g2a)
    f3a[ii] = abs(comptrap(f3,a=0,b=1,n=n[ii]) - g3a)
    f4a[ii] = abs(comptrap(f4,a=0,b=2*π,n=n[ii]) - g4a)

end

pplt.close(); fig,axs = pplt.subplots(ncols=2,axwidth=2)

lgd = Dict("frame"=>false,"ncol"=>1)
axs[1].plot(n,f1a,label=L"\int_2^{10} 3x \>dx",legend="l",legend_kw=lgd)
axs[1].plot(n,f2a,label=L"\int_0^\pi \sin(x) \>dx",legend="l")
axs[1].plot(n,f3a,label=L"\int_0^1 e^{2\cos(2\pi x)} \>dx",legend="l")
axs[1].plot(n,f4a,label=L"\int_0^{2\pi} |\cos(x)| \>dx",legend="l")
axs[1].plot([1,10000],10. .^[-1,-9],label=L"\alpha = 2",legend="l",lw=1,c="k",linestyle=":")
axs[1].format(urtitle="(a)")

nvec = vcat(1:100,100,200,500,1000,2000,5000,10000)
ln = length(nvec)
σe = zeros(ln)

for ni = 1 : ln

    A,ψ,xc,xt,σ = ExteriorMatrices1(nvec[ni])
    σa = A \ ψ
    σe[ni] = sum(abs.(σa-σ)) / nvec[ni]

end

axs[2].plot(nvec,σe,s=5,label=L"-\frac{1}{\pi}(\frac{I}{2}+f(\theta))",legend="r",legend_kw=lgd)
axs[2].scatter(nvec,σe,s=5,c="k")
axs[2].format(urtitle="(b)")

for ax in axs
    ax.format(xscale="log",yscale="log",xformatter="log",yformatter="log")
    ax.format(
        ylim=10. .^(-20,1),xlabel="n",ylabel=L"|I - $\hat{I}$|",
        suptitle="Error Decay"
    )
end

fig.savefig(plotsdir("prj3_qn0bqn1b.png"),transparent=false,dpi=300)