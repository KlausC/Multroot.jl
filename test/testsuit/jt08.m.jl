function jt08()
#
#  test polynomial suggested by Jenkins and Traub
#
    p = poly([-1.0;-1;-1;-1;-1])
    z = [-1 5]
    p, PolyZeros(z)
end
