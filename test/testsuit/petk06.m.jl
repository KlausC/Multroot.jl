function petk06()
#
# M. Petkovic testing polynomials, page 146
#
    y = [-1.0*[1,1,1,1],3*[1,1,1],-im,-im];

    p1 = poly(y);
    p2 = [1.0,-2,5];
    p2 = conv(p2,p2);
    p = conv(p1,p2);
    z = [-1.0 4; 3 3; -im 2; 1+2*im 2; 1-2*im 2];
    
    p, PolyZeros(z)
end
