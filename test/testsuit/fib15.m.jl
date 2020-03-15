function fib15()
#
#  test polynomial suggested by Goedecker
#
    n = 15;
    p = fib(n)
    z = [-.9115849521765970-.1896359497733640*im;
         -.9115849521765970+.1896359497733640*im;
         -.7620785395636261-.5388639564255547*im;
         -.7620785395636261+.5388639564255547*im;
         -.4858038058618722-.8023902123527991*im;
         -.4858038058618722+.8023902123527991*im;
         -.1248130033927766-.9371029958793902*im;
         -.1248130033927766+.9371029958793902*im;
         .2657642744875577-.9183694775911604*im;
         .2657642744875577+.9183694775911604*im;
         .6253568882929779-.7424114494678755*im;
         .6253568882929779+.7424114494678755*im;
         .8931744004970845-.4240640341420445*im;
         .8931744004970845+.4240640341420445*im;
         1.999969475434503];

    z = [z ones(n)];
    p, PolyZeros(z)
end
