function [ fd2 ] = func_fd2( n, x, beta )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fd2 = x.^n .* exp( -beta*x.^2 ) ./ ( 1.0 + exp( x ) );

end

