function igyp02(m)
#
#  generalization of Igarash and Ypma
# 
   p = poly([10*(1+i)*ones(m); -1*ones(10-m)])
   z = [10.0*(1+i) m; -1 10-m];
   p, PolyZeros(z)
end
