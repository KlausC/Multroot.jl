function bt03()
#
# Brugnano and Trigiante
#
    p = poly([im*ones(5);-im*ones(5);0.5im*ones(4);-0.5im*ones(4);0.75im;-0.75im])
    z = [im 5; -im  5; 0.5im  4; -0.5im  4; .75im  1; -0.75im  1];
    
    p, PolyZeros(z)
end
