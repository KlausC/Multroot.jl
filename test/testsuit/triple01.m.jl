function triple01()
#
#  test polynomial suggested by Goedecker
#
    p = triple(5;5;5)
    z = [[0.9;1;1.1]  [5;5;5]]

    p, PolyZeros(z)
end
