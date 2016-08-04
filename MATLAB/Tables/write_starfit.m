function [  ] = write_starfit( zz, mx, dmx, desc, fname )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

el = build_element_symbol();
mx(mx<=-99)= 0.0;


fid = fopen( fname, 'w' );

fprintf( fid, '10002\n' );
fprintf( fid, '%s\n', desc{1} );
fprintf( fid, 'Harris et al. (2016)\n' );
fprintf( fid, 'CCSN yields\n' );
fprintf( fid, '5\n' );
fprintf( fid, '\n' );
fprintf( fid, '%d\n', length(zz(mx~=0.0)) );

fform = '%-5s %+6.2f %+6.2f\n';
for inuc = 1:length(zz)
    if ( mx(inuc) ~= 0.0 )
        nname = sprintf( '%-5s', strrep(el(zz(inuc)+1,:),' ','') );
        fprintf( fid, fform, nname, mx(inuc), max(dmx(inuc),0.01) );
    end
end
fprintf( fid, '-' );
fclose(fid);

end