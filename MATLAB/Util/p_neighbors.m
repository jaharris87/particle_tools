function [ plist ] = p_neighbors( p_id, p_grid )
% Return the particle ID's corresponding to the initial immediate neighbors of particle 'p_id'

ncols = size(p_grid,1);
nrows = size(p_grid,2);

mycol = mod(p_id-1,ncols)+1;
myrow = ceil(p_id/ncols);

mycolm1 = max( mycol-1, 1 );
mycolp1 = min( mycol+1, ncols );

myrowm1 = max( myrow-1, 1 );
myrowp1 = min( myrow+1, nrows );

plist = unique( [ p_grid(mycol,myrowm1), p_grid(mycolm1:mycolp1,myrow)', p_grid(mycol,myrowp1) ] );

end