function jt06()
#
#  test polynomial suggested by Jenkins and Traub
#
    p = poly([.1,1.001,.998,1.00002,.99999])
    z = [.1;1.001;.998;1.00002;.99999];
    z = [z ones(5)]
    p, PolyZeros(z)
end
