function jt04()
#
#  test polynomial suggested by Jenkins and Traub
#
   p = poly([0.1;0.1;0.1;0.5;0.6;0.7]);
   z = [0.5  1; 0.6  1; 0.7  1; 0.1  3];
   p, PolyZeros(z)
end
