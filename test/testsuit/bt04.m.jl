function bt04()
#
# Brugnano and Trigiante
#
    p = poly([1;1;1; -1*[1;1;1;1]; (.5+im)*[1;1;1]; (.5-im)*[1;1;1];
                      .5*(1+im)*[1;1]; .5*(1-im)*[1;1]])

    z = [1 3; -1 4; .5+im 3; .5-im 3; .5*(1+im) 2; .5*(1-im) 2]

    p, PolyZeros(z)
end
