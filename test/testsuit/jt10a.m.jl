function jt10a()
#
#  test polynomial suggested by Jenkins and Traub
#
    a = 1e3;
    p = poly([a,1,1/a])
    z = [[a,1,1/a] ones(3)]
    p, PolyZeros(z)
end
