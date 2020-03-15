function jt03()
#
#  test polynomial suggested by Jenkins and Traub
#
    z = 1 ./ exp10.(1:8)
    p = poly(z)
    z = [z ones(8)]
    p, PolyZeros(z)
end
