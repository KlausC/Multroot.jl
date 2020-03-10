function inex(digits)
#
#  test polynomial suggested by Goedecker
#
    p = poly([(10/11)*[1;1;1;1;1];(20/11)*[1;1;1];(30/11)*[1;1]])
    p = round.(10^digits*p) ./ 10^digits
    z = reverse([10/11 5; 20/11 3; 30/11 2])
    p, PolyZeros(z)
end
