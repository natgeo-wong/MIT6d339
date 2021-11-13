using StatsBase

f1(x) = (6*x + 2*x^2) * exp(x)
u1(x) = 2 * x * (1-x) * exp(x)

function f2(x)
    if x <= 0.5
          return 2
    else; return 0
    end
end

function u2(x)
    if x <= 0.5
          return x * (3 - 4*x) / 4
    else; return (1-x) / 4
    end
end

f3(x) = 3*x
u3(x) = - x^3/2 + 0.5*x

p1norm(u2,u1,nx) = sum(abs.(u2.-u1)) / nx
p2norm(u2,u1) = rmsd(u2,u1)
p∞norm(u2,u1) = maximum(abs.(u2.-u1))