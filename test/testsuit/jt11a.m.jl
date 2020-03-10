function jt11a()
#
#  test polynomial suggested by Jenkins and Traub
#
    m = 15;
    y = exp.((im*pi) .* (1-m:m-1) ./ (2*m))
    z = 0.9*exp((im*pi) .*(m:3*m) ./ (2*m))
    p1 = poly(y)
    p2 = poly(z)
    p = conv(p1,p2)
    z = [[y; z] ones(4*m)]
    p, PolyZeros(z)
end
