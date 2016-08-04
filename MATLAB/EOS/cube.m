function y = cube( aa, bb )

y = ( sq( aa, bb ) + aa ).^(1/3) - ( max( sq( aa, bb ) - aa, 0.0 ) ).^(1/3);

end