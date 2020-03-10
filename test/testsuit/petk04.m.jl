function petk04()
#
# M. Petkovic testing polynomials, p139
#
    z = [3.0, -1*[1,1,1], 2*im*[1,1,1], (-2+im)*[1,1], (-2-im)*[1,1],
         (2+im)*[1,1], (2-im)*[1,1]];
    p = poly(z);
    z = [3 1; -1 3;  2*im 3; -2+im 2; -2-im 2; 2+im 2; 2-im 2];

    p, PolyZeros(z)
end
