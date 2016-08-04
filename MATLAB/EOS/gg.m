function y = gg(xx)

y = 8.0 * xx.^3 .* ( sqrt( xx.^2 + 1.0 ) - 1.0 ) - f(xx);

end