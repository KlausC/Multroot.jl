function uhlig02()
#
#  test polynomial used by F. Uhlig
#
     a = 0.001;
     p1 = [1;0;0;0;-a^4]
     p2 = reverse(poly([a;a;a;a]).a)
     p = conv(p1,p2)
     z = [ -0.00000000000000 + 0.00100000000000im
         -0.00000000000000 - 0.00100000000000im
         -0.00100000000000 + 0.00000000000000im
         0.00100000000000 + 0.00000000000000im]
    z = [z [1;1;1;5]]

    p, PolyZeros(z)
end
