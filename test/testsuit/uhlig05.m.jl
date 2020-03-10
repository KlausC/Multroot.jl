function uhlig05()
#
#  test polynomial used by F. Uhlig
#
    p = poly([-1*ones(6); 2; 2]);
    z = [-1.0 6; 2 2];

    p, PolyZeros(z)
end
