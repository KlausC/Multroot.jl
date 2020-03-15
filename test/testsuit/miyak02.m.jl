function miyak02()
#
#  test polynomial suggested by Goedecker
#  square of Fibocacci polynomial
#
    p = poly([(1.1+1.1im)*ones(4); (3.2+2.3im)*[1;1]; 2.1+1.5im])
    p = reverse(conv(p.a, p.a))
    z = [[(1.1+1.1im);(3.2+2.3im);2.1+1.5im] 2*[4;2;1]]
    p, PolyZeros(z)
end
