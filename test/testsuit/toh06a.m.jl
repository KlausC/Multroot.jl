function toh06a()
#
#  generalization of K.C. Toh and L. N. Trefethen
# 
   p = ones(21)
   z = [0.95557280578614 + 0.29475517441090im
        0.95557280578614 - 0.29475517441090im
        0.82623877431600 + 0.56332005806362im
        0.82623877431600 - 0.56332005806362im
        0.62348980185873 + 0.78183148246803im
        0.62348980185873 - 0.78183148246803im
        0.36534102436640 + 0.93087374864420im
        0.36534102436640 - 0.93087374864420im
        0.07473009358642 + 0.99720379718118im
        0.07473009358642 - 0.99720379718118im
        -0.98883082622513 + 0.14904226617617im
        -0.98883082622513 - 0.14904226617617im
        -0.90096886790242 + 0.43388373911756im
        -0.90096886790242 - 0.43388373911756im
        -0.73305187182983 + 0.68017273777092im
        -0.73305187182983 - 0.68017273777092im
        -0.50000000000000 + 0.86602540378444im
        -0.50000000000000 - 0.86602540378444im
        -0.22252093395631 + 0.97492791218182im
        -0.22252093395631 - 0.97492791218182im];
    z = [z ones(20)]

    p, PolyZeros(z)
end
