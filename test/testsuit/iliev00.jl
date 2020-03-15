function iliev(k)
#
#  generalization of Iliev example:
#
#    (x-1)^k (x-2)^2k (x-3)^4k
#  
    p = poly([ones(k);2ones(2k);3ones(3k)])
    z = [1.0 k; 2 2k; 3 4k]
    p, PolyZeros(z)
end
