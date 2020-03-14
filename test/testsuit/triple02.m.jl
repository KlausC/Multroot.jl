function triple02()
#
#  test polynomial suggested by Goedecker
#
    p = triple(10, 10, 10)
    z = [[0.9;1;1.1]  [10; 10; 10]]

    p, PolyZeros(z)
end
