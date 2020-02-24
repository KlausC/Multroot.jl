function bt01()
#
# Brugnano and Trigiante
#
    p = poly([ones(1,6),-1,-1,i,i,i,-i,-i,-i,2]);
    z = [1,6; -1, 2; i, 3; -i, 3; 2, 1];
    
    if norm(imag(z(:,1))) == 0 
        @printf("                 roots         multiplicities\n");
        @printf("\n");
        @printf("%25.15f \t \t \t %3g \n", 0.0, z');
    else
        @printf("                 roots ")
        @printf("   \t\t\t\t\t\t     multiplicities\n");
        @printf("\n");
        @printf("%22.15f + %22.15f i \t \t %3g \n", 0.0, 0.0, [real(z(:,1)),imag(z(:,1)),z(:,2)]');
    end;    
   p,z
 end