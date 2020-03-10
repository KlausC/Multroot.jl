function toh06c()
#
#  generalization of K.C. Toh and L. N. Trefethen
# 
    p = ones(6)
    p = conv(p,p)
    p = conv(p,p);
    z = [0.50000000000000 + 0.86602540378444im
         0.50000000000000 - 0.86602540378444im
        -1.00000000000000                    
        -0.50000000000000 + 0.86602540378444im
        -0.50000000000000 - 0.86602540378444im]
    z = [z 4*ones(5)]

    p, PolyZeros(z)
end
