function [ dFD2_dbeta ] = func_dFD2_dbeta( n, x, beta )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

dFD2_dbeta = - x.^n .* exp( -beta*x.^2 ) ./ ( 1.0 + exp( x ) ).^2;

end

