function y = f(xx)

y = xx .* ( 2.0 * xx.^2 - 3.0 ) .* sqrt( xx.^2 - 1.0 ) + 3.0 * log( xx + sqrt( 1.0 + xx.^2 ) );

end