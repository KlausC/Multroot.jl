function petk03()
#
# M. Petkovic testing polynomials, p134
#
    p1 = [1.0;-1-99*im/70];
    p2 = [1;2;3];
    p2 = conv(p2,p2);
    p = conv(p1,p2);
    p = conv(p,[1;1]);
    z = [-1.00000000000000 + 1.41421356237309im 2;
         -1.00000000000000 - 1.41421356237309im 2;
         -1-99*im/70 1;
         -1 1];
    p, PolyZeros(z)
end
