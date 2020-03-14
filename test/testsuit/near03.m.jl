function near03()
#
#  test polynomial suggested by Z. Zeng
#
    e = 0.001;
    p = poly([(1-e)*ones(20);ones(20);-0.5*[1;1;1;1;1]]);
    z = [1-e 20; 1 20; -0.5 5];
    p, PolyZeros(z)
end
