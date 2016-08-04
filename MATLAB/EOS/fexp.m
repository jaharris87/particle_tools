function y = fexp(xx)

expmax = 300.0;

y = exp( min( expmax, max( -expmax, xx ) ) );

end