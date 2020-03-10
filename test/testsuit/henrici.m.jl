function henrici()
#
# M. Petkovic testing polynomials, page 123
#
    p = poly([-4.1;-3.8;-2.05;-1.85;1.95;2.15;3.9;4.05])
    z = [[-4.1;-3.8;-2.05;-1.85;1.95;2.15;3.9;4.05] ones(8)];
    
    p, PolyZeros(z)
end
