function toh01()
#
#  generalization of K.C. Toh and L. N. Trefethen
# 
    k = [(-10:9)...]
    z = 2.0 .* (k .+ 0.5) ./ 19 + im .* sin.((2*pi) .* (k .+ 0.5) ./ 19)
    p = poly(z)
    z = [z ones(20)]

    p, PolyZeros(z)
end
