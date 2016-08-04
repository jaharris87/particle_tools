function [ rad, theta, vr, rho, temp, ye, flxtot, nu_temp ] = extra_particles( nrows, ncols, rad_init_all, time, grid_edge, grid_vr, grid_rho, grid_temp, grid_ye, p_rad, p_flxtot, p_nu_temp )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

nextra   = nrows * ncols;
nrowsold = length(rad_init_all) - nrows;
rad_init = rad_init_all(nrowsold+1:end);
nnu      = size(p_flxtot,2);
 
% rad  = zeros( length(time), length(nextra) );
% vr   = zeros( length(time), length(nextra) );
% rho  = zeros( length(time), length(nextra) );
% temp = zeros( length(time), length(nextra) );

% Calculate zone centers
grid_center = 0.5 * ( grid_edge(1:end-1) + grid_edge(2:end) );

% Get initial quantities of particle via zone-centered interpolation
theta_init = acos( ( ncols - 2*(1:ncols) + 1 ) / ncols );

vr_init     = interp1( grid_center, grid_vr, rad_init );
rho_init    = interp1( grid_center, grid_rho, rad_init );
temp_init   = interp1( grid_center, grid_temp, rad_init );
ye_init     = interp1( grid_center, grid_ye, rad_init );

flxtot_init  = repmat( p_flxtot, 1, 1, nrows );
nu_temp_init = repmat( p_nu_temp, 1, 1, nrows );
    
% Assume particles fall in with constant velocity
vr_row  = repmat( vr_init, length(time), 1 );
vr      = reshape( repmat( vr_row, ncols, 1 ), length(time), [] );

rad_row = repmat( rad_init, length(time), 1 ) + vr_row .* repmat( time, 1, nrows );
rad     = reshape( repmat( rad_row, ncols, 1 ), length(time), [] );

theta_row = repmat( theta_init, length(time), 1 );
theta     = repmat( theta_row, 1, nrows );

% Evolve density and temperature assuming homologous expansion
rho_row  = repmat( rho_init, length(time), 1 ) .* ( repmat( rad_init, length(time), 1 ) ./ rad_row ) .^ 3;
rho      = reshape( repmat( rho_row, ncols, 1 ), length(time), [] );

temp_row = repmat( temp_init, length(time), 1 ) .* ( repmat( rad_init, length(time), 1 ) ./ rad_row );
temp     = reshape( repmat( temp_row, ncols, 1 ), length(time), [] ); 
temp     = temp * 1e-9;

% Assume Ye doesn't change
ye_row = repmat( ye_init, length(time), 1 );
ye     = reshape( repmat( ye_row, ncols, 1 ), length(time), [] );

flxtot  = zeros( length(time), nnu, nextra );
nu_temp = zeros( length(time), nnu, nextra );

for i = 1:nnu
    % Scale neutrino fluxes with radius
    flxtot_row    = flxtot_init(:,i,:) .* reshape( ( repmat( p_rad, 1, nrows ) ./ rad_row ) .^ 2, length(time), 1, nrows );
    flxtot(:,i,:) = reshape( repmat( flxtot_row, ncols, 1, 1 ), length(time), 1, [] );

    % Assume neutrino temperatures do not change
    nu_temp_row    = nu_temp_init(:,i,:);
    nu_temp(:,i,:) = reshape( repmat( nu_temp_row, ncols, 1, 1 ), length(time), 1, [] );
end

end

