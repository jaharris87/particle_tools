function [  ] = write_inab_xnet( nname, y, desc, fname )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

ny = length(y);
y(y<=1.0e-99)= 0.0;

fid = fopen( fname, 'w' );

fprintf( fid, '%s\n', desc );
fform = strcat( repmat('% 5s% 23.15E ',1,4),'\n' );
output = cell(2,ny);
for i = 1:ny
    output{1,i} = strrep(nname(i,:),' ','');
    output{2,i} = y(i);
end
fprintf( fid, fform, output{:,:} );
fclose(fid);

end