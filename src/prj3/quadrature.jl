function comptrap(
    f :: Function;
    a :: Real,
    b :: Real,
    n :: Int
)

    g = (f(a) + f(b))/2

    for ii = 1 : (n-1)

        g += f(a + ii * (b-a)/n)

    end

    return g * (b-a) / n
    
end