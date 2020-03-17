function bt02()
#
# Brugnano and Trigiante
#
    p = poly([ones(10);-1;-1;im;-im;2])
    z = [1.0 10; -1 2; im 1; -im 1; 2 1];
    
    p, PolyZeros(z)
end
