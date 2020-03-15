function fl(k)
#
#  generalization of Farmer-Loizou example:
#
#    (x-1)^4k * (x-2)^3k * (x-3)^2K * (x-4)^k
#
    p = poly([ones(4*k);2*ones(3*k);3*ones(2*k);4*ones(k)])
    z = [1.0 4k; 2 4k; 3 2k; 4 k]
    p, PolyZeros(z)
end
