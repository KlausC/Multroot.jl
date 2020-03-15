function petk01()
#
# M. Petkovic testing polynomials (p. 109)
#
    p1 = poly([-1.0; -1; 3; 3; 3; -im*[1;1;1;1]]);
    p1 = reverse(p1.a)
    p2 = [1;-2;5];
    p2 = conv(p2,p2);
    p = conv(p1,p2);
    z = [-1 2; 3  3; -im  4; 1+2*im  2; 1-2*im  2];
    p, PolyZeros(z)
end
