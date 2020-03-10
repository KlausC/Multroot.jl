function jt07a()
#
#  test polynomial suggested by Jenkins and Traub
#
   a = 10.0^(-10);
   p = poly([.001; .01; .1; .1+a*i; .1-a*i; 1; -10])
   z = [.001; .01; .1; .1+a*im; .1-a*im; 1; -10]
   z = [z ones(7)]
   y = [-10  1; 1  1; 0.01  1; 0.001  1; 0.1  3];
   p, PolyZeros(z), PolyZeros(y)
end
