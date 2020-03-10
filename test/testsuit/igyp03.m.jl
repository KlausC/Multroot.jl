function igyp03()
#
#  generalization of Igarash and Ypma
# 
    p = poly([10.0*(1+im)*[1;1;1];1;im;2;2im;3;4im;5])
    z = [10*(1+im) 3; 1 1; im 1; 2 1; 2im 1; 3 1; 4im 1; 5 1];
    p, PolyZeros(z)
end
