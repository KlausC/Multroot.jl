function twin01()

   k = 4
   p = poly([-0.2*ones(k); 0.39*ones(k); 0.40*ones(k)])
   z = [-0.2 k; 0.39 k; 0.4 k];

   p, PolyZeros(z)
end
