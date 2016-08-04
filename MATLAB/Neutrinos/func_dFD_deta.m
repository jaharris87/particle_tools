function [ dFD_deta ] = func_dFD_deta( n, x, eta )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

dFD_deta = x.^n .* exp( x - eta ) ./ ( 1.0 + exp( x - eta ) ).^2;

end

