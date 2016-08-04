%% pp_mass_itime
% Interpolate isotope masses to specified time for each particle
function [ varargout ] = pp_mass_interp( time, plist, p_data, pp_time, pp_xn, tracer_mass )
%% Input Arguments
%   Variable        Type        Dimension   Units   Description
%   --------        -------     ---------   -----   ------------------------------------------------------------
%   time            Float        M          s       1D array of elapsed time
%     >=0
%
%   plist           Integer      N                  1D array containing particle IDs to be used
%     >0                                                (eg. 1:4000 or 1:40:4000 or [70,543,800,...])
%
%   p_data          Integer     (M, N)              Particle fate flags at each timestep in 'time' from find_particle_fates routine
%     >=0
% 
%   pp_time         Cell         N          s       Contains the elapsed times for each particle ts file
%     >=0
% 
%   pp_xn           Cell         N                  Contains time history of mass fraction for every isotope for each particle ts file
%     >=0
%
%   tracer_mass     Float        N          g       Mass of each particle
%     >0
%
%% Output Arguments
%   Variable            Type        Dimension   Units   Description
%   ---------           -------     ---------   -----   ------------------------------------------------------------
%   pp_mx               Cell         N          M_sun   Mass of each isotope at for each particle at every timestep in 'pp_time'
%
%   pp_mxinterp_posvr   Integer     (M, K)      M_sun   Time history of total mass of each isotope in unbound particles with positive radial velocity
%
%   pp_mxinterp_unb     Integer     (M, K)      M_sun   Time history of total mass of each isotope in unbound particles
%
%% Optional Input Arguments
%   Flag            Variable        Type        Dimension   Units   Description                                                     Default
%   ---------       --------        -------     ---------   -----   ------------------------------------------------------------    ---------------------------
%
%% Initialization

load_constants;

% Create an instance of the inputParser class.
p = inputParser;
p.FunctionName = 'pp_mass_interp';
p.KeepUnmatched = true;

% Define required inputs
p.addRequired('time', ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'vector', ...
                                         'real', ...
                                         '>=', 0.0}));
p.addRequired('plist', ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'vector', ...
                                         'integer', ...
                                         '>', 0}));
p.addRequired('p_data', ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'2d', ...
                                         'integer', ...
                                         '>=', 0}));
p.addRequired('pp_time', ...
              @(x)validateattributes(x, {'cell'}, ...
                                        {'vector'}));
p.addRequired('pp_xn', ...
              @(x)validateattributes(x, {'cell'}, ...
                                        {'vector'}));
p.addRequired('tracer_mass', ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'vector', ...
                                         'real', ...
                                         '>', 0.0}));

% Parse, validate, and assign required input arguments.
p.parse( time, plist, p_data, pp_time, pp_xn, tracer_mass );
time            = p.Results.time;
plist           = p.Results.plist;
p_data          = p.Results.p_data;
pp_time         = p.Results.pp_time;
pp_xn           = p.Results.pp_xn;
tracer_mass     = p.Results.tracer_mass;

% Get the dimensions from input
nparticles = length(plist);
nisotopes  = length( pp_xn{1}(:,1) );

% Create the output arrays
pp_mx             = cell(1,nparticles);
pp_mxinterp_posvr = zeros(length(time),nisotopes);
pp_mxinterp_unb   = zeros(length(time),nisotopes);
    
% Temporary array for interpolating
pp_mx_tmp         = zeros(length(time),nparticles);

%% Calculate isotope masses and mass total for unbound and positive radial velocity masks

for i = 1:nparticles

    pp_mx{i} = pp_xn{i} .* repmat( tracer_mass(i), nisotopes, length( pp_time{i} ) ) / M_solar;

%     posvr_mask = p_data(:,i)==3;
%     unb_mask   = p_data(:,i)==3 | p_data(:,i)==2;
    
    % Do the interpolation only if enough data points
%     if length( pp_time{i} ) > 1
%         for j = 1:nisotopes
%             pp_mx_tmp(:,i)                  = interp1( pp_time{i}, pp_mx{i}(j,:), time, 'linear', pp_mx{i}(j,end) );
%             pp_mxinterp_posvr(posvr_mask,j) = pp_mxinterp_posvr(posvr_mask,j) + pp_mx_tmp(posvr_mask,i);
%             pp_mxinterp_unb(unb_mask,j)     = pp_mxinterp_unb(unb_mask,j) + pp_mx_tmp(unb_mask,i);
%         end
%     else
%         pp_mx_tmp2 = repmat(pp_mx{i},1,length(time))';
%         pp_mxinterp_posvr(posvr_mask,:) = pp_mxinterp_posvr(posvr_mask,:) + pp_mx_tmp2(posvr_mask,:);
%         pp_mxinterp_unb(unb_mask,:) = pp_mxinterp_unb(unb_mask,:) + pp_mx_tmp2(unb_mask,:);
%     end
    
end

%% Finalize

varargout{1} = pp_mx;
if nargout > 1
    varargout{2} = pp_mxinterp_posvr;
    if nargout > 2
        varargout{3} = pp_mxinterp_unb;
    end
end

end

