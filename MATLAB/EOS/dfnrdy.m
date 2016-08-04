function y = dfnrdy(yy)

y = yy.^1.5 * ( 1.0 + pi^2 / ( 8.0 * yy.^2 ) + 7.0 * pi^4 / (640.0 * yy.^4 ) );

end