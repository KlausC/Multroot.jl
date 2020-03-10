function bt01()
#
# Brugnano and Trigiante
#
    p = poly([ones(6);-1;-1;im;im;im;-im;-im;-im;2])
    z = [1 6; -1  2; im  3; -im  3; 2  1];
    p, PolyZeros(z)
 end
