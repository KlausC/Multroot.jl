function petk02()
#
# M. Petkovic testing polynomials, p118
#
    p = poly([1;1;-im;-im;-im;5*im;5*im;-5*im;-5*im]);
    z = [1 2; -im 3; 5*im 2; -5*im 2];
    
    p, PolyZeros(z)
end
