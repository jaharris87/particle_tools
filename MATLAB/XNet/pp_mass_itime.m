%% pp_mass_itime
% Interpolate isotope masses to specified time for each particle
function [ varargout ] = pp_mass_itime( time, time_interp, plist, pp_time, pp_mx )
%% Input Arguments
%   Variable        Type        Dimension   Units   Description
%   --------        -------     ---------   -----   ------------------------------------------------------------
%   time            Float        M          s       1D array of elapsed time
%     >=0
% 
%   time_interp     Float        1          s       Time at which to interpolate
%     >=0
%
%   plist           Integer      N                  1D array containing particle IDs to be used
%     >0                                                (eg. 1:4000 or 1:40:4000 or [70,543,800,...])
% 
%   pp_time         Cell         N          s       Contains the elapsed times for each particle ts file
%     >=0
% 
%   pp_mx           Cell         N          M_sun   Contains time history of mass for every isotope for each particle ts file
%     >=0
%
%% Output Arguments
%   Variable        Type        Dimension   Units   Description
%   ---------       -------     ---------   -----   ------------------------------------------------------------
%   pp_mx_itime     Float       (K, N)      M_sun   Mass of each isotope at for each particle at time(itime)
%
%   itime           Integer      1                  Index of 'time' array closest to 'time_interp'
%
%% Optional Input Arguments
%   Flag            Variable        Type        Dimension   Units   Description                                                     Default
%   ---------       --------        -------     ---------   -----   ------------------------------------------------------------    ---------------------------
%
%% Initialization

load_constants;

% Create an instance of the inputParser class.
p = inputParser;
p.FunctionName = 'pp_mass_itime';
p.KeepUnmatched = true;

% Define required inputs
p.addRequired('time', ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'vector', ...
                                         'real', ...
                                         '>=', 0.0}));
p.addRequired('time_interp', ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'real', ...
                                         '>=', 0.0}));
p.addRequired('plist', ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'vector', ...
                                         'integer', ...
                                         '>', 0}));
p.addRequired('pp_time', ...
              @(x)validateattributes(x, {'cell'}, ...
                                        {'vector'}));
p.addRequired('pp_mx', ...
              @(x)validateattributes(x, {'cell'}, ...
                                        {'vector'}));

% Parse, validate, and assign required input arguments.
p.parse( time, time_interp, plist, pp_time, pp_mx );
time            = p.Results.time;
time_interp     = p.Results.time_interp;
plist           = p.Results.plist;
pp_time         = p.Results.pp_time;
pp_mx           = p.Results.pp_mx;
                                 
% Get the dimensions from input
nparticles = length(plist);
nisotopes  = length( pp_mx{1}(:,1) );

% Initialize to zero
pp_mx_itime = zeros( nisotopes, nparticles );

%% Interpolate masses to time closest to 'time_interp' in 'time' array

% Find the index in 'time' closest to 'time_interp'
[ ~, itime ] = min( abs( time - time_interp ) );

% Do the interpolation only if enough data points
for i = 1:nparticles
    if length(pp_time{i}) > 2
        for j = 1:nisotopes
            pp_mx_itime(j,i) = interp1( pp_time{i}, pp_mx{i}(j,:), time(itime), 'linear', pp_mx{i}(j,end) );
        end
    else
        pp_mx_itime(:,i) = pp_mx{i}(:,end);
    end
end

%% Finalize

varargout{1} = pp_mx_itime;
if nargout > 1
    varargout{2} = itime;
end

end

