function triple04()
#
#  test polynomial suggested by Goedecker
#
    p = triple(25, 15, 10)
    z = [[0.9;1;1.1]  [25; 15; 10]]

    p, PolyZeros(z)
end
