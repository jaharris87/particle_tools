function [  ] = write_table_z( zz, mx, desc, fname )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

el = build_element_symbol();
mx(mx<=1.0e-99)= 0.0;

nmx = size(mx,2);

fid = fopen( fname, 'w' );

% fform = '% 33s\n';
% fprintf( fid, fform, 'Mass (solar)' );
fform = strcat( '% 5s% 4s',repmat('% 24s',1,nmx),'\n' );
fprintf( fid, fform, 'name', 'Z', desc{:} );
fform = strcat( '% 5s% 4d',repmat('% 24.15E',1,nmx),'\n' );
for inuc = 1:length(zz)
    if ( mx(inuc,1) < Inf )
        nname = sprintf( '% 5s', lower(strrep(el(zz(inuc)+1,:),' ','')) );
        fprintf( fid, fform, nname, zz(inuc), mx(inuc,1:end) );
    end
end
fclose(fid);

end