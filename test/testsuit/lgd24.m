function [p,z] = lgd24
%
%  test polynomial suggested by Goedecker
%   Legendre polynomial of degree 24
%
    p = lgd(24);
    z = [-.9951872199970214, ... 
        -.9747285559713095, ... 
        -.9382745520027328, ... 
        -.8864155270044010, ... 
        -.8200019859739029, ... 
        -.7401241915785544, ... 
        -.6480936519369756, ... 
        -.5454214713888395, ... 
        -.4337935076260451, ... 
        -.3150426796961634, ... 
        -.1911188674736163, ... 
        -.6405689286260563e-1, ... 
        .6405689286260563e-1, ... 
        .1911188674736163, ... 
        .3150426796961634, ... 
        .4337935076260451, ... 
        .5454214713888395, ... 
        .6480936519369756, ... 
        .7401241915785544, ... 
        .8200019859739029, ... 
        .8864155270044010, ... 
        .9382745520027328, ... 
        .9747285559713095, ... 
        .9951872199970214];
    z = [z',ones(24,1)];
    fprintf('Legendre polynomial, roots are sensitive \n');
        fprintf('                 roots         multiplicities\n');
        fprintf('\n');
        fprintf('%25.14f \t \t \t %3g \n', z');
      
    