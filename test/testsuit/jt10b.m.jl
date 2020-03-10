function jt10b()
#
#  test polynomial suggested by Jenkins and Traub
#
    a = 1e6;
    p = poly([a,1,1/a])
    z = [[a,1,1/a] ones(3)]
    p, PolyZeros(z)
end
