%% write_th_file
% write_th_file writes various thermodynamic profiles for all particles in a model and in parallel if proper flag is set
function [ varargout ] = get_th_tau( plist, fate_fname_base, varargin )
%% Input Arguments
%   Variable        Type        Dimension   Units   Description
%   --------        -------     ---------   -----   ------------------------------------------------------------
%   plist           Integer      N                  1D array containing particle IDs to be used
%     >0                                                (eg. 1:4000 or 1:40:4000 or [70,543,800,...])
%
%   fate_fname_base String       1                  Base filename for raw ASCII "fate" data
%
%% Output Arguments
%   Variable        Type        Dimension   Units   Description
%   ---------       -------     ---------   -----   ------------------------------------------------------------
%   tau             Float        N          s       Expansion timescale for each particle in 'plist' [s]
%
%% Optional Input Arguments
%   Flag            Variable        Type        Dimension   Units   Description                                                     Default
%   ---------       --------        -------     ---------   -----   ------------------------------------------------------------    ---------------------------
%   ChimeraNSE      chimera_flag    Logical      1                  Flag for using Chimera NSE data to set t_nse                    false
%
%   ChimeraTime     t_chim          Float        K          s       Times corresponding to Chimera NSE data                         []
%     >=0.0
%
%   Parallel        parallel_flag   Logical      1                  Flag for writing profiles in parallel                           false
%
%   RadiusNSE       radius_nse      Float       (K, L)      cm      Radii for NSE transition in Chimera simulation for each L       []
%     >=0                                                               angular zones at every K timestep
%
%   TempNSE         temp_nse        Float        1          GK      NSE transition temperature from which to begin profile          8.0
%     >=0                                                               ( if >0.0, transition temperature determined by
%                                                                         find_nse_temp if chimera_flag == false or
%                                                                         find_nse_chimera if chimera_flag == true )
%
%   ThetaNSE        theta_nse       Float        L          rad     Angular position of each Chimera radial ray                     []
%     >=0
%
%   TimeStop        t_stop_array    Float        1          s       Final time for profiles and time from which to perform          time(end)
%     >=0                                                               extrapolation if extrap_flag == true
%
%   AdiabaticTol    adiabatic_tol   Float        1                  Tolerance for "constant" adiabatic criteria
%     >0                                                                ( T * rho^{-1/3} = constant )
%                                                                       ( Default: 0.05 )
%
%   MinExtrapTime   time_extrap_min Float        1          s       Minimum duration of being within adiabatic_tol
%     >=0                                                               required to extrapolate
%                                                                       ( Default: 0.0025 )
%
%   MaxExtrapTime   time_extrap_max Float        1          s       Maximum lookback time for determination of
%     >=0                                                               extrapolation parameters
%                                                                       ( Default: 0.15 )
%
%   SGWinTau        sgwin_tau       Float        1          s       Window size of Savitzky-Golay smoothing of
%     >0                                                                expansion timescale
%                                                                       ( Default: 0.025 )
%
%   SGOrderTau      sgorder_tau     Integer      1                  Order of Savitzky-Golay smoothing of
%     >=1                                                               expansion timescale
%                                                                       ( Default: 6 )
%
%   SGWinAdi        sgwin_adi       Float        1          s       Window size of Savitzky-Golay smoothing of
%     >0                                                                "constant" adiabatic criteria
%                                                                       ( Default: 0.05 )
%
%   SGOrderAdi      sgorder_adi     Integer      1                  Order of Savitzky-Golay smoothing of
%     >=1                                                               "constant" adiabatic criteria
%                                                                       ( Default: 2 )
%
%   MinPosTau       min_jumpdiff    Integer      1                  Minimum number of consecutive positive values
%     >=1                                                               of expansion timescale to be used in calculation
%                                                                       ( Default: 5 )
%
%   TauTol          tau_rho_tol     Float        1                  Criteria for eliminating large variations in
%     >0                                                                expansion timescale
%
%   TauMin          tau_rho_min     Float        1          s       Minimum value allowed for expansion timescale
%     >=0                                                               to be used in calculation
%                                                                       ( Default: 0.01 )
%
%   TauMax          tau_rho_max     Float        1          s       Maximum value allowed for expansion timescale
%     >=0                                                               to be used in calculation
%                                                                       ( Default: 1.0 )
%
%% Initialization

% Create an instance of the inputParser class.
p = inputParser;
p.FunctionName = 'get_th_tau';
p.KeepUnmatched = true;

% Define required inputs
p.addRequired('plist', ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'vector', ...
                                         'integer', ...
                                         '>', 0}));
p.addRequired('fate_fname_base', ...
              @(x)validateattributes(x, {'char'}, ...
                                        {'row'}));
% Define optional inputs
p.addOptional('ChimeraNSE', false, ...
              @(x)validateattributes(x, {'logical'}, ...
                                        {'scalar'}));
p.addOptional('ChimeraTime', [], ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'vector', ...
                                         'real', ...
                                         '>=', 0.0}));
p.addOptional('Parallel', false, ...
              @(x)validateattributes(x, {'logical'}, ...
                                        {'scalar'}));
p.addOptional('RadiusNSE', [], ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'2d', ...
                                         'real', ...
                                         '>=', 0.0}));
p.addOptional('TempNSE', 8.0, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'real', ...
                                         '>=', 0.0}));
p.addOptional('ThetaNSE', [], ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'2d', ...
                                         'real', ...
                                         '>=', 0.0}));
p.addOptional('TimeStop', 0.0, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'real', ...
                                         '>=', 0.0}));
                                     
p.addOptional('AdiabaticTol', 0.05, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'real', ...
                                         '>', 0.0}));
p.addOptional('MinExtrapTime', 0.0025, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'real', ...
                                         '>=', 0.0}));
p.addOptional('MaxExtrapTime', 0.15, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'real', ...
                                         '>=', 0.0}));
p.addOptional('SGWinTau', 0.025, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'real', ...
                                         '>', 0.0}));
p.addOptional('SGOrderTau', 6, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'integer', ...
                                         '>=', 1}));
p.addOptional('SGWinAdi', 0.05, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'real', ...
                                         '>', 0.0}));
p.addOptional('SGOrderAdi', 2, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'integer', ...
                                         '>=', 1}));
p.addOptional('MinPosTau', 5, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'integer', ...
                                         '>=', 1}));
p.addOptional('TauMin', 0.01, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'real', ...
                                         '>=', 0.0}));
p.addOptional('TauMax', 1.0, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'real', ...
                                         '>=', 0.0}));
                                     
% Parse, validate, and assign input arguments.
p.parse( plist, fate_fname_base, varargin{:} );
plist           = p.Results.plist;
fate_fname_base = p.Results.fate_fname_base;
chimera_flag    = p.Results.ChimeraNSE;
t_chim          = p.Results.ChimeraTime;
parallel_flag   = p.Results.Parallel;
radius_nse      = p.Results.RadiusNSE;
temp_nse        = p.Results.TempNSE;
theta_nse       = p.Results.ThetaNSE;
t_stop_array    = p.Results.TimeStop;
adiabatic_tol   = p.Results.AdiabaticTol;
time_extrap_min = p.Results.MinExtrapTime;
time_extrap_max = p.Results.MaxExtrapTime;
sgwin_tau       = p.Results.SGWinTau;
sgorder_tau     = p.Results.SGOrderTau;
sgwin_adi       = p.Results.SGWinAdi;
sgorder_adi     = p.Results.SGOrderAdi;
min_jumpdiff    = p.Results.MinPosTau;
tau_rho_min     = p.Results.TauMin;
tau_rho_max     = p.Results.TauMax;

% Override temp_nse if chimera_flag is set to true
if chimera_flag
    temp_nse = 0.0;
end

plist_all = plist;

if size( t_stop_array, 1 ) == length(plist_all)
    t_stop_array = t_stop_array';
end

if size( t_stop_array, 2 ) == 1
    t_stop_array(:,plist_all) = t_stop_array(:,1);
end

% Initialize output variables
tau = zeros( size(t_stop_array) );

% Set flags for recording output variables
if nargout > 0
    keep_tau = true;
else
    keep_tau = false;
end
    
%% Write Files

if parallel_flag
    
    % Launch a parallel job
    poolobj = parpool;

    % Initialize progress monitor
%   N = length(plist_all);
%   percentage_update = 0.05;
%   do_debug = 0;
%   run_javaaddpath = 1;
%   show_execution_time = 1;
%   ppm = ParforProgressStarter2('write_th_file', N, percentage_update, do_debug, run_javaaddpath, show_execution_time);
    
    parfor id = 1:length(plist_all)
        
        % Get particle ID and output it to screen
        p_id = plist_all(id);
        
        % Load particle data from files
        [ time, radius, theta, ~, ~, temp, density, ~, ~, ~, ~, ~, ~, ~, ~ ] = read_tracer_fate( fate_fname_base, p_id );
        
        t_stop = t_stop_array(:,id);
        t_stop(t_stop == 0.0) = time(end);
        
        time = time( time <= max(t_stop) );
        radius = radius( time <= max(t_stop) );
        density = density( time <= max(t_stop) );
        temp = temp( time <= max(t_stop) );
        
        % Check to see if profile should begin at a constant NSE transition temperature or calculated
        if temp_nse > 0.0

            % Extrapolation to temp_final (if temp_final > 0.0), transitioning from NSE at temperature temp_nse
            ptau = find_expansion_timescale( time, temp, density, p_id, ...
                'TimeStop', t_stop, ...
                'TempNSE', temp_nse, ...
                'AdiabaticTol', adiabatic_tol, ...
                'MinExtrapTime', time_extrap_min, ...
                'MaxExtrapTime', time_extrap_max, ...
                'SGWinTau', sgwin_tau, ...
                'SGOrderTau', sgorder_tau, ...
                'SGWinAdi', sgwin_adi, ...
                'SGOrderAdi', sgorder_adi, ...
                'MinPosTau', min_jumpdiff, ...
                'TauMin', tau_rho_min, ...
                'TauMax', tau_rho_max );
            
        else
            
            % Replicate CHIMERA NSE behavior
            if chimera_flag
                
                % Determine time of NSE transition by comparing particle data to Chimera NSE data
                [ t_nse, ~, ~ ] = find_nse_chimera( time, radius, theta, t_chim, radius_nse, theta_nse, p_id );
                
            else
                
                % Determine time of NSE transition as a function of density and temperature
                [ t_nse, ~, ~ ] = find_nse_temp( time, density, temp, p_id );
                
            end
            
            % Extrapolation to temp_final (if temp_final > 0.0), transitioning from NSE at time t_nse
            ptau = find_expansion_timescale( time, temp, density, p_id, ...
                'TimeStart', t_nse, ...
                'TimeStop', t_stop, ...
                'AdiabaticTol', adiabatic_tol, ...
                'MinExtrapTime', time_extrap_min, ...
                'MaxExtrapTime', time_extrap_max, ...
                'SGWinTau', sgwin_tau, ...
                'SGOrderTau', sgorder_tau, ...
                'SGWinAdi', sgwin_adi, ...
                'SGOrderAdi', sgorder_adi, ...
                'MinPosTau', min_jumpdiff, ...
                'TauMin', tau_rho_min, ...
                'TauMax', tau_rho_max );
            
        end
        
        % Keep timescale and profile data if needed
        if keep_tau
            tau(:,id) = ptau;
        end
        
        % Update progress monitor
%       ppm.increment(id);
        
    end
    
    % Clean-up for progress monitor
%   delete(ppm);
    
    % Delete the parallel job
    delete(poolobj);
    
    
else
    
    for id = 1:length(plist_all)
        
        % Get particle ID and output it to screen
        p_id = plist_all(id);
        disp(p_id);
        
        % Load particle data from files
        [ time, radius, theta, ~, ~, temp, density, ~, ~, ~, ~, ~, ~, ~, ~ ] = read_tracer_fate( fate_fname_base, p_id );
        
        t_stop = t_stop_array(:,id);
        t_stop(t_stop == 0.0) = time(end);
        
        time = time( time <= max(t_stop) );
        radius = radius( time <= max(t_stop) );
        density = density( time <= max(t_stop) );
        temp = temp( time <= max(t_stop) );
        
        % Check to see if profile should begin at a constant NSE transition temperature or be calculated
        if temp_nse > 0.0

            % Extrapolation to temp_final (if temp_final > 0.0), transitioning from NSE at temperature temp_nse
            ptau = find_expansion_timescale( time, temp, density, p_id, ...
                'TimeStop', t_stop, ...
                'TempNSE', temp_nse, ...
                'AdiabaticTol', adiabatic_tol, ...
                'MinExtrapTime', time_extrap_min, ...
                'MaxExtrapTime', time_extrap_max, ...
                'SGWinTau', sgwin_tau, ...
                'SGOrderTau', sgorder_tau, ...
                'SGWinAdi', sgwin_adi, ...
                'SGOrderAdi', sgorder_adi, ...
                'MinPosTau', min_jumpdiff, ...
                'TauMin', tau_rho_min, ...
                'TauMax', tau_rho_max );
            
        else
            
            % Replicate CHIMERA NSE behavior
            if chimera_flag
                
                % Determine time of NSE transition by comparing particle data to Chimera NSE data
                [ t_nse, ~, ~ ] = find_nse_chimera( time, radius, theta, t_chim, radius_nse, theta_nse, p_id );
                
            else
                
                % Determine time of NSE transition as a function of density and temperature
                [ t_nse, ~, ~ ] = find_nse_temp( time, density, temp, p_id );
                
            end
            
            % Extrapolation to temp_final (if temp_final > 0.0), transitioning from NSE at time t_nse
            ptau = find_expansion_timescale( time, temp, density, p_id, ...
                'TimeStart', t_nse, ...
                'TimeStop', t_stop, ...
                'AdiabaticTol', adiabatic_tol, ...
                'MinExtrapTime', time_extrap_min, ...
                'MaxExtrapTime', time_extrap_max, ...
                'SGWinTau', sgwin_tau, ...
                'SGOrderTau', sgorder_tau, ...
                'SGWinAdi', sgwin_adi, ...
                'SGOrderAdi', sgorder_adi, ...
                'MinPosTau', min_jumpdiff, ...
                'TauMin', tau_rho_min, ...
                'TauMax', tau_rho_max );
            
        end
        
        % Keep timescale and profile data if needed
        if keep_tau
            tau(:,id) = ptau;
        end
        
    end
    
end

%% Finalize

varargout{1} = tau;

end