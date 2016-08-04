function [ dfbeta_dbeta ] = dfbeta_dbeta( beta )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

ul = 50;

g2 = integral(@(x)func_fd2(2,x,beta),0,ul);
g3 = integral(@(x)func_fd2(3,x,beta),0,ul);
g4 = integral(@(x)func_fd2(4,x,beta),0,ul);

dg2 = integral(@(x)func_dFD2_dbeta(2,x,beta),0,ul);
dg3 = integral(@(x)func_dFD2_dbeta(3,x,beta),0,ul);
dg4 = integral(@(x)func_dFD2_dbeta(4,x,beta),0,ul);

dfbeta_dbeta = (dg2*g4 + g2*dg4 - 2*g2*g4*dg3/g3)/g3^2;

end

