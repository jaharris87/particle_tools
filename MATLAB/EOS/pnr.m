function y = pnr(yy)

y = yy.^2.5 * ( 1.0 + 5.0 * pi^2 / ( 8.0 * yy.^2 ) - 7.0 * pi^4 / ( 384.0 * yy.^4 ) );

end