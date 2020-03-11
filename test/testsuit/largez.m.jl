function largez(k::Integer)
#
#  test polynomial suggested by Z Zeng
#
    p = [ 
         1.0
        -0.7
        -0.19
         0.177
        -0.7364
        -0.43780
        -0.952494
        -0.2998258
        -0.00322203
        -0.328903811
        -0.4959527435
        -0.9616679762
         0.4410459281
         0.1090273141
         0.6868094008
         0.0391923826
         0.0302248540
         0.6603775863
        -0.1425784968
        -0.3437618593
         0.4357949015
        ];
    z = [
         0.5 + im
         0.5 - im
        -1 + 0.2im
        -1 - 0.2im
        -0.1 + im
        -0.1 - im
         0.8 + 0.6im
         0.8 - 0.6im
        -0.7 + 0.7im
        -0.7 - 0.7im
         1.4
        -0.4 + 0.9im
        -0.4 - 0.9im
         0.9
        -0.8 + 0.3im
        -0.8 - 0.3im
         0.3 + 0.8im
         0.3 - 0.8im
         0.6 + 0.4im
         0.6 - 0.4im
        ]
    
    mult = 1
    for i = 1:k
        p = conv(p, p)
        mult *= 2
    end

    z = [z mult*ones(20)]
    p, PolyZeros(z)
end

largez01() = largez(0)
largez02() = largez(1)
largez03() = largez(2)
largez04() = largez(3)
largez05() = largez(4)
largez06() = largez(5)

