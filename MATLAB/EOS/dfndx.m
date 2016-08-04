function y = dfndx( xx, beta )

y = 3.0 * ffn( xx, beta ) ./ xx + (pi/beta)^2 * ( 2.0 - 4.0 * ( 2.0 * xx.^2 + 1.0 ) / xx.^2 );

end