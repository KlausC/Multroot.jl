function toh05()
#
#  generalization of K.C. Toh and L. N. Trefethen
# 
    k = [(-19:0)...]
    z = exp2.(k)
    p = poly(z)
    z = [z ones(20)]

    p, PolyZeros(z)
end
