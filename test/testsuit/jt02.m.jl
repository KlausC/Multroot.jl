function jt02()
#
#  test polynomial suggested by Jenkins and Traub
#
    p = poly([1.0:17])
    z = [1.0:17]
    z = [z ones(17)]
    p, PolyZeros(z)
end
