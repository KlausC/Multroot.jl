function toh03()
#
#  generalization of K.C. Toh and L. N. Trefethen
# 
    k = [(1:20)...]
    z = 10/11 .- exp2.(-k)
    p = poly(z)
    z = [z ones(20)]

    p, PolyZeros(z)
end
