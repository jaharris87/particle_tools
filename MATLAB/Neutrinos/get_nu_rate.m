function [ lamba_nu ] = get_nu_rate( flxtot, lapse )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Load Physical Constants
load_constants;

unu = unui ./ lapse;

imax = size(lapse,1);

for i = 1:imax
    phi_eff = flxtot .* ( unu + dmnp_mev ) .^2 ./ ( unu .^ 2 );
end

end

