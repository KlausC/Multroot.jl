function petk07()
#
# M. Petkovic testing polynomials, page 147
#
    y = [1.0*[1,1,1],-2+im, -2-im, 5*im,5*im, -5*im, -5*im];
    p = poly(y);
    z = [1.0 3; -2+im 1;  -2-im 1;  5*im 2; -5*im 2];
    
    p, PolyZeros(z)
end
