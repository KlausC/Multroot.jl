function jt07d()
#
#  test polynomial suggested by Jenkins and Traub
#
    a = 1e-7
    z = [.001; .01; .1; .1+a*im; .1-a*im; 1; -10];
    p = poly(z)
    z = [z ones(7)]
    y = [-10  1; 1  1; 0.01  1; 0.001  1; 0.1  3];

    p, PolyZeros(z), PolyZeros(y)
end
