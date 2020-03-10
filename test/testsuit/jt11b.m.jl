function jt11b()
#
#  test polynomial suggested by Jenkins and Traub
#
    m = 20;
    y = exp.((im*pi) .* (1-m:m-1) ./ (2*m))
    z = 0.9*exp.((im*pi) .* (m:3*m) ./ (2*m))
    u = [y; z]
    p = poly(u)
    z = [[y; z] ones(4*m)]
    p, PolyZeros(z)
end
