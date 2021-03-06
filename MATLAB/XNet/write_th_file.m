%% write_th_file
% write_th_file writes various thermodynamic profiles for all particles in a model and in parallel if proper flag is set
function [ varargout ] = write_th_file( plist, time_bounce, fate_fname_base, temp_fname_base, th_fname_base, varargin )
%% Input Arguments
%   Variable        Type        Dimension   Units   Description
%   --------        -------     ---------   -----   ------------------------------------------------------------
%   plist           Integer      N                  1D array containing particle IDs to be used
%     >0                                                (eg. 1:4000 or 1:40:4000 or [70,543,800,...])
%
%   time_bounce     Float        1          s       Time of core-bounce
%     >=0
%
%   fate_fname_base String       1                  Base filename for raw ASCII "fate" data
%
%   temp_fname_base String       1                  Base filename for raw ASCII "temp" data
%
%   th_fname_base   String       1                  Base filename for ASCII thermo profiles
%
%% Output Arguments
%   Variable        Type        Dimension   Units   Description
%   ---------       -------     ---------   -----   ------------------------------------------------------------
%   tau             Float        N          s       Expansion timescale for each particle in 'plist' [s]
%
%   th_out          Cell         N                  Thermodynamic profile data for each particle in 'plist'
%
%% Optional Input Arguments
%   Flag            Variable        Type        Dimension   Units   Description                                                     Default
%   ---------       --------        -------     ---------   -----   ------------------------------------------------------------    ---------------------------
%   ChangeMin       change_min      Float        1                  Min change in thermo conditions for writing data to profile     0.001
%     >=0
%
%   ChimeraNSE      chimera_flag    Logical      1                  Flag for using Chimera NSE data to set t_nse                    false
%
%   ChimeraTime     t_chim          Float        K          s       Times corresponding to Chimera NSE data                         []
%     >=0.0
%
%   ExtraRows       nrows_extra     Integer      M                  Number of additional particle rows for which to generate        0
%     >=0                                                               homologous contraction profiles
%
%   Extrapolate     extrap_flag     Logical      1                  Flag for extrapolating profile to temp_final                    false
%
%   ExtrapTemp      temp_final      Float        1          GK      Final temperature in extrapolation                              0.5
%     >0
%
%   ExtrapTime      time_final      Float        1          Final time in extrapolation [s]
%     >0                                                        ( Default: time to reach temp_final )
%
%   GridFile        grid_fname      String       1                  Filename for WH07 initial progenitor file
%
%   NCols           ncols           Integer      1                  Total original number of particle columns for model             40
%     >0
%
%   NRows           nrows           Integer      1                  Total original number of particle rows for model                ceil( N / ncols )
%     >=0
%
%   OutputTimestep  output_tdel     Float        1          s       Output to profile via constant timestep if >0.0                 0.0
%     >=0
%
%   Parallel        parallel_flag   Logical      1                  Flag for writing profiles in parallel                           false
%
%   RadiusNSE       radius_nse      Float       (K, L)      cm      Radii for NSE transition in Chimera simulation for each L       []
%     >=0                                                               angular zones at every K timestep
%
%   RowRadii        prad_rows       Float       (nrows + M) cm      Initial radii for all particle rows, including extra rows       []
%     >=0                                                               (must be defined if extra rows are specified)
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
%   NuTimeStop      nu_time_stop    Float        1          s       Time at which to zero neutrino-fluxes
%     =>0                                                               ( Default: 10.0 )
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
p.FunctionName = 'write_th_file';
p.KeepUnmatched = true;

% Define required inputs
p.addRequired('plist', ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'vector', ...
                                         'integer', ...
                                         '>', 0}));
p.addRequired('time_bounce', ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'real', ...
                                         '>=', 0.0}));
p.addRequired('fate_fname_base', ...
              @(x)validateattributes(x, {'char'}, ...
                                        {'row'}));
p.addRequired('temp_fname_base', ...
              @(x)validateattributes(x, {'char'}, ...
                                        {'row'}));
p.addRequired('th_fname_base', ...
              @(x)validateattributes(x, {'char'}, ...
                                        {'row'}));
% Define optional inputs
p.addOptional('ChangeMin', 0.001, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'real', ...
                                         '>=', 0.0}));
p.addOptional('ChimeraNSE', false, ...
              @(x)validateattributes(x, {'logical'}, ...
                                        {'scalar'}));
p.addOptional('ChimeraTime', [], ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'vector', ...
                                         'real', ...
                                         '>=', 0.0}));
p.addOptional('ExtraRows', 0, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'integer', ...
                                         '>=', 0}));
p.addOptional('Extrapolate', false, ...
              @(x)validateattributes(x, {'logical'}, ...
                                        {'scalar'}));
p.addOptional('ExtrapTemp', 0.5, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'real', ...
                                         '>', 0.0}));
p.addOptional('ExtrapTime', 0.0, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'real', ...
                                         '>', 0.0}));
p.addOptional('GridFile', '', ...
              @(x)validateattributes(x, {'char'}, ...
                                        {'row'}));
p.addOptional('NCols', 40, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'integer', ...
                                         '>', 0}));
p.addOptional('NRows', 0, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'integer', ...
                                         '>=', 0}));
p.addOptional('OutputTimestep', 0.0, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'real', ...
                                         '>=', 0.0}));
p.addOptional('Parallel', false, ...
              @(x)validateattributes(x, {'logical'}, ...
                                        {'scalar'}));
p.addOptional('RadiusNSE', [], ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'2d', ...
                                         'real', ...
                                         '>=', 0.0}));
p.addOptional('RowRadii', [], ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'vector', ...
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
p.addOptional('NuTimeStop', 10.0, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'real', ...
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
p.parse( plist, time_bounce, fate_fname_base, temp_fname_base, th_fname_base, varargin{:} );
plist           = p.Results.plist;
time_bounce     = p.Results.time_bounce;
fate_fname_base = p.Results.fate_fname_base;
temp_fname_base = p.Results.temp_fname_base;
th_fname_base   = p.Results.th_fname_base;
change_min      = p.Results.ChangeMin;
chimera_flag    = p.Results.ChimeraNSE;
t_chim          = p.Results.ChimeraTime;
nrows_extra     = p.Results.ExtraRows;
extrap_flag     = p.Results.Extrapolate;
temp_final      = p.Results.ExtrapTemp;
time_final      = p.Results.ExtrapTime;
grid_fname      = p.Results.GridFile;
ncols           = p.Results.NCols;
nrows           = p.Results.NRows;
output_tdel     = p.Results.OutputTimestep;
parallel_flag   = p.Results.Parallel;
radius_nse      = p.Results.RadiusNSE;
prad_rows       = p.Results.RowRadii;
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
nu_time_stop    = p.Results.NuTimeStop;

% Calculate number of rows if value is not supplied, assuming plist contains all of the particles in the model
if nrows == 0
    nrows = ceil( length(plist) / ncols );
end

% Override temp_nse if chimera_flag is set to true
if chimera_flag
    temp_nse = 0.0;
end

% Deal with any additional particles
if nrows_extra > 0
    
    % Calculate the number of extra particles and assign particle IDs
    nrows_all = nrows + nrows_extra;
    n_extra   = nrows_extra * ncols;
    p_extra   = (nrows * ncols) + (1:n_extra);
    plist_all = [ plist, p_extra ];
    
    % Get ID and neutrino data for particle in middle of last row of particles which will be scaled appropriately for the extra particles
    p_copy    = nrows * ncols - floor(ncols/2) - 1;
    [ t_copy, rad_copy, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~ ] = read_tracer_fate( fate_fname_base, p_copy );
    [ ~, ~, ~, ~, flxtot_copy, nu_temp_copy ] = read_tracer_temp( temp_fname_base, p_copy );
    
    % Read the initial grid data which will be interpolated to give the initial state for the extra particles
    [ ~, rad_grid, vr_grid, rho_grid, temp_grid, ~, ~, ~, ~, ~, ye_grid ] = import_grid_wh07( grid_fname );
    rad_grid = [ 0.0; rad_grid ];
    
    % Calculate trajectories of extra particles under the assumption of homologous contraction
    t_extra = t_copy;
    [ rad_extra, theta_extra, vr_extra, rho_extra, temp_extra, ye_extra, flxtot_extra, nu_temp_extra ] ...
        = extra_particles( nrows_all, ncols, prad_rows, t_copy, rad_grid, vr_grid, rho_grid, temp_grid, ye_grid, rad_copy, flxtot_copy, nu_temp_copy );
    
else
    
    t_extra = [];
    rad_extra = [];
    theta_extra = [];
    vr_extra = [];
    rho_extra = [];
    temp_extra = [];
    ye_extra = [];
    flxtot_extra = [];
    nu_temp_extra = [];
    
    plist_all = plist;
    
end

if length(t_stop_array) == 1 && length(plist_all) > 1
    t_stop_array(plist_all) = t_stop_array(1);
end

% Initialize output variables
tau = zeros( size(plist_all) );
th_out = cell( size(plist_all) );

% Set flags for recording output variables
if nargout > 0
    keep_tau = true;
    if nargout > 1
        keep_th = true;
    else
        keep_th = false;
    end
else
    keep_tau = false;
    keep_th  = false;
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
%         disp(p_id);
        
        % Check if this particle was originally in the model
        if p_id <= nrows * ncols
            
            % Load particle data from files
            [ time, radius, theta, v_rad, ~, temp, density, ye, ~, ~, ~, ~, ~, ~, ~ ] = read_tracer_fate( fate_fname_base, p_id );
            [ ~, ~, ~, ~, flxtot, nu_temp ] = read_tracer_temp( temp_fname_base, p_id );
            
        else
            
            % Load homologous contraction data
            time    = t_extra;
            radius  = rad_extra(:,p_id);
            theta   = theta_extra(:,p_id);
            v_rad   = vr_extra(:,p_id);
            density = rho_extra(:,p_id);
            temp    = temp_extra(:,p_id);
            ye      = ye_extra(:,p_id);
            flxtot  = flxtot_extra(:,:,p_id);
            nu_temp = nu_temp_extra(:,:,p_id);
            
        end
        
        if t_stop_array(id) == 0.0
            t_stop = time(end);
        else
            t_stop = t_stop_array(id);
        end
        
        time = time( time <= t_stop );
        radius = radius( time <= t_stop );
        v_rad = v_rad( time <= t_stop );
        density = density( time <= t_stop );
        temp = temp( time <= t_stop );
        ye = ye( time <= t_stop );
        flxtot = flxtot( time <= t_stop, : );
        nu_temp = nu_temp( time <= t_stop, : );
        
        % Check to see if profile should begin at a constant NSE transition temperature or calculated
        if temp_nse > 0.0

            % Extrapolation to temp_final (if temp_final > 0.0), transitioning from NSE at temperature temp_nse
            [ ptau, pth_out ] = plot_expansion( time, temp, density, ye, flxtot, nu_temp, radius, v_rad, p_id, time_bounce, ...
                'Extrapolate', extrap_flag, ...
                'ExtrapTemp', temp_final, ...
                'ExtrapTime', time_final, ...
                'TempNSE', temp_nse, ...
                'OutputTimestep', output_tdel, ...
                'TimeStop', t_stop, ...
                'ChangeMin', change_min, ...
                'AdiabaticTol', adiabatic_tol, ...
                'MinExtrapTime', time_extrap_min, ...
                'MaxExtrapTime', time_extrap_max, ...
                'NuTimeStop', nu_time_stop, ...
                'SGWinTau', sgwin_tau, ...
                'SGOrderTau', sgorder_tau, ...
                'SGWinAdi', sgwin_adi, ...
                'SGOrderAdi', sgorder_adi, ...
                'MinPosTau', min_jumpdiff, ...
                'TauMin', tau_rho_min, ...
                'TauMax', tau_rho_max, ...
                'PlotProfile', false, ...
                'SavePlot', false, ...
                'WriteProfile', true, ...
                'ProfileBaseName', th_fname_base );
            
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
            [ ptau, pth_out ] = plot_expansion( time, temp, density, ye, flxtot, nu_temp, radius, v_rad, p_id, time_bounce, ...
                'Extrapolate', extrap_flag, ...
                'ExtrapTemp', temp_final, ...
                'ExtrapTime', time_final, ...
                'TimeStart', t_nse, ...
                'TimeStop', t_stop, ...
                'OutputTimestep', output_tdel, ...
                'ChangeMin', change_min, ...
                'AdiabaticTol', adiabatic_tol, ...
                'MinExtrapTime', time_extrap_min, ...
                'MaxExtrapTime', time_extrap_max, ...
                'NuTimeStop', nu_time_stop, ...
                'SGWinTau', sgwin_tau, ...
                'SGOrderTau', sgorder_tau, ...
                'SGWinAdi', sgwin_adi, ...
                'SGOrderAdi', sgorder_adi, ...
                'MinPosTau', min_jumpdiff, ...
                'TauMin', tau_rho_min, ...
                'TauMax', tau_rho_max, ...
                'PlotProfile', false, ...
                'SavePlot', false, ...
                'WriteProfile', true, ...
                'ProfileBaseName', th_fname_base );
            
        end
        
        % Keep timescale and profile data if needed
        if keep_tau
            tau(id) = ptau;
        end
        if keep_th
            th_out{id} = pth_out;
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
        
        % Check if this particle was originally in the model
        if p_id <= nrows * ncols
            
            % Load particle data from files
            [ time, radius, theta, v_rad, ~, temp, density, ye, ~, ~, ~, ~, ~, ~, ~ ] = read_tracer_fate( fate_fname_base, p_id );
            [ ~, ~, ~, ~, flxtot, nu_temp ] = read_tracer_temp( temp_fname_base, p_id );
            
        else
            
            % Load homologous contraction data
            time    = t_extra;
            radius  = rad_extra(:,p_id);
            theta   = theta_extra(:,p_id);
            v_rad   = vr_extra(:,p_id);
            density = rho_extra(:,p_id);
            temp    = temp_extra(:,p_id);
            ye      = ye_extra(:,p_id);
            flxtot  = flxtot_extra(:,:,p_id);
            nu_temp = nu_temp_extra(:,:,p_id);
            
        end
        
        if t_stop_array(id) == 0.0
            t_stop = time(end);
        else
            t_stop = t_stop_array(id);
        end
        
        time = time( time <= t_stop );
        radius = radius( time <= t_stop );
        v_rad = v_rad( time <= t_stop );
        density = density( time <= t_stop );
        temp = temp( time <= t_stop );
        ye = ye( time <= t_stop );
        flxtot = flxtot( time <= t_stop, : );
        nu_temp = nu_temp( time <= t_stop, : );
        
        % Check to see if profile should begin at a constant NSE transition temperature or be calculated
        if temp_nse > 0.0

            % Extrapolation to temp_final (if temp_final > 0.0), transitioning from NSE at temperature temp_nse
            [ ptau, pth_out ] = plot_expansion( time, temp, density, ye, flxtot, nu_temp, radius, v_rad, p_id, time_bounce, ...
                'Extrapolate', extrap_flag, ...
                'ExtrapTemp', temp_final, ...
                'ExtrapTime', time_final, ...
                'TempNSE', temp_nse, ...
                'OutputTimestep', output_tdel, ...
                'TimeStop', t_stop, ...
                'ChangeMin', change_min, ...
                'AdiabaticTol', adiabatic_tol, ...
                'MinExtrapTime', time_extrap_min, ...
                'MaxExtrapTime', time_extrap_max, ...
                'NuTimeStop', nu_time_stop, ...
                'SGWinTau', sgwin_tau, ...
                'SGOrderTau', sgorder_tau, ...
                'SGWinAdi', sgwin_adi, ...
                'SGOrderAdi', sgorder_adi, ...
                'MinPosTau', min_jumpdiff, ...
                'TauMin', tau_rho_min, ...
                'TauMax', tau_rho_max, ...
                'PlotProfile', false, ...
                'SavePlot', false, ...
                'WriteProfile', true, ...
                'ProfileBaseName', th_fname_base );
            
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
            [ ptau, pth_out ] = plot_expansion( time, temp, density, ye, flxtot, nu_temp, radius, v_rad, p_id, time_bounce, ...
                'Extrapolate', extrap_flag, ...
                'ExtrapTemp', temp_final, ...
                'ExtrapTime', time_final, ...
                'TimeStart', t_nse, ...
                'TimeStop', t_stop, ...
                'OutputTimestep', output_tdel, ...
                'ChangeMin', change_min, ...
                'AdiabaticTol', adiabatic_tol, ...
                'MinExtrapTime', time_extrap_min, ...
                'MaxExtrapTime', time_extrap_max, ...
                'NuTimeStop', nu_time_stop, ...
                'SGWinTau', sgwin_tau, ...
                'SGOrderTau', sgorder_tau, ...
                'SGWinAdi', sgwin_adi, ...
                'SGOrderAdi', sgorder_adi, ...
                'MinPosTau', min_jumpdiff, ...
                'TauMin', tau_rho_min, ...
                'TauMax', tau_rho_max, ...
                'PlotProfile', false, ...
                'SavePlot', false, ...
                'WriteProfile', true, ...
                'ProfileBaseName', th_fname_base );
            
        end
        
        % Keep timescale and profile data if needed
        if keep_tau
            tau(id) = ptau;
        end
        if keep_th
            th_out{id} = pth_out;
        end
        
    end
    
end

%% Finalize

varargout{1} = tau;
if nargout > 1
    varargout{2} = th_out;
end

end