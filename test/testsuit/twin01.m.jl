function twin01()

   k = 4;
   p = poly([-0.2*ones(1,k),0.39*ones(1,k),0.40*ones(1,k)])
   z = [-0.2, k; 0.39, k; 0.4, k];
   @printf("                 roots         multiplicities\n");
        @printf("\n");
        @printf("%25.15f \t \t \t %3g \n", z');p,z
end