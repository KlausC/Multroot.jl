function jt05()
#
#  test polynomial suggested by Jenkins and Traub
#
    p = poly([0.1*ones(4);0.2*ones(3);0.3;0.3;0.4])
    z = [0.4  1; 0.3  2; 0.2  3; 0.1  4];
    p, PolyZeros(z)
end
