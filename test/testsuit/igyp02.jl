function igyp02(m)
#
#  generalization of Igarash and Ypma
# 
   p = poly([10*(1+im)*ones(m); -1*ones(10-m)])
   z = [10.0*(1+im) m; -1 10-m];
   p, PolyZeros(z)
end
