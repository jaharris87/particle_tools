function [ dfeta_deta ] = dfeta_deta( eta )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

f2 = integral(@(x)func_fd(2,x,eta),0,inf);
f3 = integral(@(x)func_fd(3,x,eta),0,inf);
f4 = integral(@(x)func_fd(4,x,eta),0,inf);

df2 = integral(@(x)func_dFD_deta(2,x,eta),0,inf);
df3 = integral(@(x)func_dFD_deta(3,x,eta),0,inf);
df4 = integral(@(x)func_dFD_deta(4,x,eta),0,inf);

dfeta_deta = (df2*f4 + f2*df4 + f2*f4 - 2*f2*f4*df3/f3)/f3^2;

end

