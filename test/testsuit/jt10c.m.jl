function jt10c()
#
#  test polynomial suggested by Jenkins and Traub
#
    a = 1e9;
    p = poly([a,1,1/a])
    z = [[a,1,1/a] ones(3)]
    p, PolyZeros(z)
end
