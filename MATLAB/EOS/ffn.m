function y = ffn( xx, beta )

y = xx.^3 .* ( 1.0 + (pi/beta)^2 * ( 2.0 * xx.^2 + 1.0 ) / ( 2.0 * xx.^4 ) + 6.0/40.0 * (pi/beta)^2 / xx.^8 );

end