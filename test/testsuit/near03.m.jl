function near03()
#
#  test polynomial suggested by Z. Zeng
#
    e = 0.001;
    p = poly([(1-e)*ones(1,20),ones(1,20),-0.5*[1,1,1,1,1]]);
    z = [1-e, 20; 1, 20; -0.5, 5];
    if norm(imag(z(:,1))) == 0 
        @printf("                 roots         multiplicities\n");
        @printf("\n");
        @printf("%25.15f \t \t \t %3g \n", z');
    else
        @printf("                 roots ")
        @printf("   \t\t\t\t\t\t     multiplicities\n");
        @printf("\n");
        @printf("%22.15f + %22.15f i \t \t %3g \n",             [real(z(:,1)),imag(z(:,1)),z(:,2)]');
    end;    
    p,z
end