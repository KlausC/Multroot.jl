function uhlig04()
#
#  test polynomial used by F. Uhlig
#
    r = [(3/11)*ones(12); 11/3; 11/3; (2*im/7)*ones(4); (2.5+im/4)*ones(2); 1/8]
    p = poly(r);
    z = [3/11 12; 11/3 2; 2*im/7 4; 2.5+im/4 2; 1/8 1];

    p, PolyZeros(z)
end
