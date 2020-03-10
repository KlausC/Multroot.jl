function lgd05()
#
#  test polynomial suggested by Goedecker
#   Legendre polynomial of degree 5
#
    p = lgd(5);
    z = [0,
         0.90617984593866,
         0.53846931010568,
         -0.90617984593866,
         -0.53846931010568
        ];
    z = [z ones(5)]
    p, PolyZeros(z)
end
