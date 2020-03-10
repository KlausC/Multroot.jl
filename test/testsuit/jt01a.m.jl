function jt01a()
#
#  test polynomial suggested by Jenkins and Traub
#
   a = 10.0^10
   p = [1;-1;-a^2;a^2]
   z = [a 1;-a 1;1 1]
   p, PolyZeros(z)
end
