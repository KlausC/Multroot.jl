function triple03()
#
#  test polynomial suggested by Goedecker
#
    p = triple(18, 10, 16)
    z = [[0.9;1;1.1]  [18; 10; 16]]

    p, PolyZeros(z)
end
