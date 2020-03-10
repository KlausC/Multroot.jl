function fibsq04()
#
#  test polynomial suggested by Goedecker
#
    n = 4;
    p = fib(n)
    p = conv(p,p)
    z = [-.7748041132154339;
         -.07637893113374573-.8147036471703865*im; -.07637893113374573+.8147036471703865*im;
         1.927561975482925]
    z = [z  2*ones(n)]
    p, PolyZeros(z)
end
