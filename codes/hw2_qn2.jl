using DrWatson
@quickactivate "MIT6.339"

using PyCall, LaTeXStrings
pplt = pyimport("proplot")

function uh(h::Real)

    uh = 0
    for j = 1 : (1/h)
        uh += 1 / j^2
    end

    return uh

end

eh(h::Real) = abs(uh(h)-π^2/6)

h = exp.(-10:0.1:0)

pplt.close(); fig,axs = pplt.subplots(axwidth=2)

axs[1].plot(log.(h),log.(eh.(h)))
axs[1].format(xlabel="h",ylabel=L"e_h",xlim=(-10,0),ylim=(-10,0))

fig.savefig(plotsdir("hw2_qn2.png"),transparent=false,dpi=250)