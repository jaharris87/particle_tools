function [ fd ] = func_fd( n, x, eta )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fd = x.^n ./ ( 1.0 + exp( x - eta ) );

end

