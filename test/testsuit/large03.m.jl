function large03()
#
#  test polynomial suggested by Z Zeng
#
   
    p = [ 
        1.0;
        0;
        -0.28000000000000;
        -0.59200000000000;
        -0.96890000000000;
        -0.62020000000000;
        0.54382200000000;
        0.08842680000000;
        0.43204191000000;
        -0.15503004000000;
        0.44119337520000;
        -0.65443195750000;
        -0.35708097260000;
        -0.55597791620000;
        0.63916879370000;
        0.14637576020000;
        -0.06709193010000;
        -0.03193754010000;
        0.50857405010000;
        -0.14968258190000;
        0.63273024920000 ];
    p = conv(p,p);
    p = conv(p,p);
    z = [-1.0 + 0.3im;
        -1.0 - 0.30im;
        -0.9 + 0.4im;
        -0.9 - 0.4im;
        -0.7 + 0.7im;
        -0.7 - 0.7im;
        -0.4 + 0.9im;
        -0.4 - 0.9im;
        -0.0 + 1.10000000000im;
        -0.000000000000 - 1.10000000000im;
        1.2;
        1.00000000000000;
        0.9 + 0.4im;
        0.9 - 0.4im;
        0.6 + 0.600000000000im;
        0.6 - 0.600000000000im;
        0.4 + 0.900000000000im;
        0.4 - 0.900000000000im;
        0.00000000000 + 0.8im
        0.00000000000 - 0.8im];
    z = [z 4*ones(20)]
    p, PolyZeros(z)
end
