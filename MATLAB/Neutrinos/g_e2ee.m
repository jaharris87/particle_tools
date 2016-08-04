function [ g_e2ee ] = g_e2ee( beta )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

g2 = integral(@(x)func_fd2(2,x,beta),0,inf);
g3 = integral(@(x)func_fd2(3,x,beta),0,inf);
g4 = integral(@(x)func_fd2(4,x,beta),0,inf);

g_e2ee = g4*g2/g3^2;

end