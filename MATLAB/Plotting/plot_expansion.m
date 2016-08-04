%% plot_expansion
% plot_expansion generates thermodynamic profiles for tracer particles for input into XNet
function [ varargout ] = plot_expansion( time, temp, density, ye, flxtot, nu_temp, radius, v_rad, plist, time_bounce, varargin )
%% Input Arguments
%   Variable        Type        Dimension   Description
%   --------        -------     ---------   -------------------------------------------------
%   time            Float        M          1D array of elapsed time for the data [s]
%     >=0
% 
%   temp            Float       (M, N)      2D array of temperature for all particles 
%     >0                                        in 'plist' at times in 'time' [GK]
% 
%   density         Float       (M, N)      2D array of density for all particles
%     >0                                        in 'plist' at times in 'time' [g/cm^{3}]
% 
%   ye              Float       (M, N)      2D array of electron fractions for all particles 
%     >=0                                       in 'plist' at times in 'time'
%     <=1
%
%   flxtot          Float       (M, 4, N)   3D array of neutrino number flux for all particles 
%                                               in 'plist' at times in 'time' of each neutrino
%                                               species [(# neutrinos)/cm^{3}/s]
%
%   nu_temp         Float       (M, 4, N)   3D array of neutrino temperature for all particles 
%     >0                                        in 'plist' at times in 'time' of each neutrino
%                                               species [MeV]
% 
%   radius          Float       (M, N)      2D array of radial position for all particles 
%     >=0                                       in 'plist' at times in 'time' [cm]
% 
%   v_rad           Float       (M, N)      2D array of radial velocity for all particles 
%     <c                                        in 'plist' at times in 'time' [cm/s]
% 
%   plist           Integer      N          1D array containing particle IDs to be used
%     >0                                        (eg. 1:4000 or 1:40:4000 or [70,543,800,...]
%
%   time_bounce     Float        1          Time of core-bounce [s]
%     >=0
%
%% Output Arguments
%   Variable        Type        Dimension   Description
%   ---------       -------     ---------   -------------------------------------------------
%   tau             Float        N          Expansion timescale for each particle in 'plist' [s]
%
%   th_out          Cell         N          Thermodynamic profile data for each particle in 'plist'
%
%% Optional Input Arguments
%   Flag            Variable        Type        Dimension   Description
%   ---------       --------        -------     ---------   -------------------------------------------------
%   AdiabaticTol    adiabatic_tol   Float        1          Tolerance for "constant" adiabatic criteria
%     >0                                                        ( T * rho^{-1/3} = constant )
%                                                               ( Default: 0.05 )
%
%   MinExtrapTime   time_extrap_min Float        1          Minimum duration of being within adiabatic_tol
%     >=0                                                       required to extrapolate [s]
%                                                               ( Default: 0.0025 )
%
%   MaxExtrapTime   time_extrap_max Float        1          Maximum lookback time for determination of
%     >=0                                                       extrapolation parameters [s]
%                                                               ( Default: 0.15 )
%
%   SGWinTau        sgwin_tau       Float        1          Window size of Savitzky-Golay smoothing of
%     >0                                                        expansion timescale [s]
%                                                               ( Default: 0.025 )
%
%   SGOrderTau      sgorder_tau     Integer      1          Order of Savitzky-Golay smoothing of
%     >=1                                                       expansion timescale
%                                                               ( Default: 6 )
%
%   SGWinAdi        sgwin_adi       Float        1          Window size of Savitzky-Golay smoothing of
%     >0                                                        "constant" adiabatic criteria [s]
%                                                               ( Default: 0.05 )
%
%   SGOrderAdi      sgorder_adi     Integer      1          Order of Savitzky-Golay smoothing of
%     >=1                                                       "constant" adiabatic criteria
%                                                               ( Default: 2 )
%
%   MinPosTau       min_jumpdiff    Integer      1          Minimum number of consecutive positive values
%     >=1                                                       of expansion timescale to be used in calculation
%                                                               ( Default: 5 )
%
%   TauTol          tau_rho_tol     Float        1          Criteria for eliminating large variations in
%     >0                                                        expansion timescale
%
%   TauMin          tau_rho_min     Float        1          Minimum value allowed for expansion timescale
%     >=0                                                       to be used in calculation [s]
%                                                               ( Default: 0.01 )
%
%   TauMax          tau_rho_max     Float        1          Maximum value allowed for expansion timescale
%     >=0                                                       to be used in calculation [s]
%                                                               ( Default: 1.0 )
%
%   ChangeMin       change_min      Float        1          Minimum change in thermodynamic conditions for
%     >=0                                                       writing data to profile
%                                                               ( Default: 0.001 )
%
%   PlotProfile     plot_flag       Logical      1          Flag for generating plots of temperature profile
%                                                               ( Default: true )
%
%   SavePlot        print_flag      Logical      1          Flag for saving plots of temperature profile to PDF
%                                                               ( Default: false )
%
%   WriteProfile    write_flag      Logical      1          Flag for generating ASCII thermodynamic profile
%                                                               for XNet
%                                                               ( Default: false )
%
%   PlotFolder      plot_folder     String       1          Folder in which to save plots
%                                                               ( Default: './' )
%
%   ProfileBaseName th_fname_base   String       1          Base filename for ASCII profiles on which to append
%                                                               the particle ID # to make the full filename
%                                                               ( Default: './th_profile/th-' )
%
%   Extrapolate     extrap_flag     Logical      1          Flag for extrapolating profile to temp_final
%                                                               ( Default: true )
%
%   ExtrapTemp      temp_final      Float        1          Final temperature in extrapolation [GK]
%     >0                                                        ( Default: 0.5 )
%
%   ExtrapTime      time_final      Float        1          Final time in extrapolation [s]
%     >=0                                                       ( Default: time to reach temp_final )
%
%   TimeExtend      t_extend        Float        1          Time to extend using constant final temperature if >0.0 [s]
%     =>0                                                       ( Default: 0.0 )
%
%   TempNSE         temp_nse        Float        1          Temperature from which to perform extrapolation [GK]
%     >0                                                        ( Default: 8.0 )
%
%   TempMin         temp_min        Float        1          Minimum temperature in extrapolation [GK]
%     >0                                                        ( Default: 0.02 )
%
%   TimeStart       t_start_array   Float        N or 1     Initial time for post-processing profile. [s]
%     >=0                                                       If N=1, the same value is used for all particles
%                                                               ( Default: last time when temp < temp_nse )
%
%   TimeStop        t_stop_array    Float        K          Times from which to perform extrapolation if
%     >=0                                                       extrap_flag == true [s]
%                                                               ( Default: time(end) )
%
%   OutputTimestep  output_tdel     Float        1          Output to profile via constant timestep if >0.0
%     >=0                                                       ( Default: 0.0 )
%
%   NuTimeStop      nu_time_stop    Float        1          Time at which to zero neutrino-fluxes [s]
%     =>0                                                       ( Default: 10.0 )
%                                                               
%% Initialization

c = 2.99792458E+10; % Speed of light [cm s^{-1}]

% Create an instance of the inputParser class.
p = inputParser;
p.FunctionName = 'plot_expansion';
p.KeepUnmatched = true;

% Define required inputs
p.addRequired('time', ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'vector', ...
                                         'real', ...
                                         '>=', 0.0}));
p.addRequired('temp', ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'2d', ...
                                         'real', ...
                                         '>', 0.0}));
p.addRequired('density', ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'2d', ...
                                         'real', ...
                                         '>', 0.0}));
p.addRequired('ye', ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'2d', ...
                                         'real', ...
                                         '>=', 0.0, ...
                                         '<=', 1.0}));
p.addRequired('flxtot', ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'size', [NaN,4,NaN], ...
                                         'real'}));
p.addRequired('nu_temp', ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'size', [NaN,4,NaN], ...
                                         'real', ...
                                         '>=', 0.0}));
p.addRequired('radius', ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'2d', ...
                                         'real', ...
                                         '>=', 0.0}));
p.addRequired('v_rad', ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'2d', ...
                                         'real', ...
                                         '<', c}));
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
% Define optional inputs
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
p.addOptional('ChangeMin', 0.001, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'real', ...
                                         '>=', 0.0}));
p.addOptional('ExtrapTemp', 0.5, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'real', ...
                                         '>', 0.0}));
p.addOptional('PlotProfile', true, ...
              @(x)validateattributes(x, {'logical'}, ...
                                        {'scalar'}));
p.addOptional('SavePlot', false, ...
              @(x)validateattributes(x, {'logical'}, ...
                                        {'scalar'}));
p.addOptional('WriteProfile', false, ...
              @(x)validateattributes(x, {'logical'}, ...
                                        {'scalar'}));                               
p.addOptional('PlotFolder', './', ...
              @(x)validateattributes(x, {'char'}, ...
                                        {'row'}));
p.addOptional('ProfileBaseName', './th_profile/th-', ...
              @(x)validateattributes(x, {'char'}, ...
                                        {'row'}));
p.addOptional('TempNSE', 8.0, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'real', ...
                                         '>=', 0.0}));
p.addOptional('TempMin', 0.02, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'real', ...
                                         '>=', 0.0}));
p.addOptional('TimeStart', 0.0, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'real', ...
                                         '>=', 0.0}));
p.addOptional('TimeStop', 0.0, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'real', ...
                                         '>=', 0.0}));
p.addOptional('Extrapolate', true, ...
              @(x)validateattributes(x, {'logical'}, ...
                                        {'scalar'}));
p.addOptional('ExtrapTime', 0.0, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'real', ...
                                         '>=', 0.0}));
p.addOptional('OutputTimestep', 0.0, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'real', ...
                                         '>=', 0.0}));
p.addOptional('TimeExtend', 0.0, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'real', ...
                                         '>=', 0.0}));
p.addOptional('NuTimeStop', 10.0, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'real', ...
                                         '>=', 0.0}));

% Parse, validate, and assign required input arguments.
p.parse( time, temp, density, ye, flxtot, nu_temp, radius, v_rad, plist, time_bounce, varargin{:} );
time            = p.Results.time;
temp            = p.Results.temp;
density         = p.Results.density;
ye              = p.Results.ye;
flxtot          = p.Results.flxtot;
nu_temp         = p.Results.nu_temp;
radius          = p.Results.radius;
v_rad           = p.Results.v_rad;
plist           = p.Results.plist;
time_bounce     = p.Results.time_bounce;
adiabatic_tol   = p.Results.AdiabaticTol;
time_extrap_min = p.Results.MinExtrapTime;
time_extrap_max = p.Results.MaxExtrapTime;
sgwin_tau       = p.Results.SGWinTau;
sgorder_tau     = p.Results.SGOrderTau;
sgwin_adi       = p.Results.SGWinAdi;
sgorder_adi     = p.Results.SGOrderAdi;
min_jumpdiff    = p.Results.MinPosTau;
% tau_rho_tol   = p.Results.TauTol;
tau_rho_min     = p.Results.TauMin;
tau_rho_max     = p.Results.TauMax;
change_min      = p.Results.ChangeMin;
plot_flag       = p.Results.PlotProfile;
write_flag      = p.Results.WriteProfile;
print_flag      = p.Results.SavePlot;
plot_folder     = p.Results.PlotFolder;
th_fname_base   = p.Results.ProfileBaseName;
temp_nse        = p.Results.TempNSE;
temp_min        = p.Results.TempMin;
t_start_array   = p.Results.TimeStart;
t_stop_array    = p.Results.TimeStop;
extrap_flag     = p.Results.Extrapolate;
temp_final      = p.Results.ExtrapTemp;
time_final      = p.Results.ExtrapTime;
output_tdel     = p.Results.OutputTimestep;
t_extend        = p.Results.TimeExtend;
nu_time_stop    = p.Results.NuTimeStop;

prc_min = 25; prc_max = 75;

time_extrap_max0 = time_extrap_max;

% t_stop_array = unique( t_stop_array( t_stop_array <= time(end) & t_stop_array > time_bounce + time_extrap_max ) );
if t_stop_array == 0.0
    t_stop_array = time(end);
end

if time_final > 0.0
    t_stop_array(t_stop_array > time_final) = time_final;
end

if print_flag
    plot_flag = true;
    if ~isdir(plot_folder)
        mkdir(plot_folder);
    end
end

if write_flag
    [ profile_folder, ~, ~ ] = fileparts( th_fname_base );
    if ~isdir(profile_folder)
        mkdir(profile_folder);
    end
    th_format = [ repmat( '% 15.7E', 1, 4 ), repmat( '% 12.3E', 1, 8 ), '\n' ];
end

adiabatic_string = sprintf( '%s%g%s', '${\Delta}(T{\rho}^{-1/3}) <$ ', ...
    adiabatic_tol, ...
    '$\times(T{\rho}^{-1/3})_{\mathrm{f}}$' );

if length(t_stop_array) > 6
    color_array = varycolor( length(t_stop_array) );
elseif length(t_stop_array) < 1
    color_array=[];
elseif length(t_stop_array)==1
    color_array=[0 0 1];
elseif length(t_stop_array)==2
    color_array=[0 0 1; 1 0 0];
elseif length(t_stop_array)==3
    color_array=[0 0 1; 1 0 0; 0 0.5 0];
elseif length(t_stop_array)==4
    color_array=[0 0 1; 1 0 0; 0 0 0.5; 0 0.75 0.75];
elseif length(t_stop_array)==5
    color_array=[0 0 1; 1 0 0; 0 0 0.5; 0 0.75 0.75; 0.75 0 0.75];
elseif length(t_stop_array)==6
    color_array=[0 0 1; 1 0 0; 0 0 0.5; 0 0.75 0.75; 0.75 0 0.75; 0.75 0.75 0];
end

if length(t_start_array) == 1 && length(plist) > 1
    t_start_array(plist) = t_start_array(1);
end

i_peak    = ones( size(plist) );

temp_peak = max(temp);
for id = plist
    
    if length(plist) == 1
        p_id = 1;
    else
        p_id = id;
    end
    
    if t_start_array(p_id) > 0.0 
        [~,i_nse] = min( abs( time - t_start_array(p_id) ) );
        i_peak(p_id) = i_nse;
    else
        i_nse = 1;
        if temp_peak(p_id) >= temp_nse && temp(end,p_id) < temp_nse
            i_nse = find( temp(:,p_id) > temp_nse, 1, 'last' );
        elseif temp(end,p_id) >= temp_nse
            i_nse = length(time);
        elseif temp_peak(p_id) < temp_nse
            i_nse = 1;
        end
        t_start_array(p_id) = time(i_nse);
        i_peak(p_id)        = i_nse;
    end
end
tau = zeros(length(t_stop_array),length(plist));
th_out = 0;

j = 1;
for id = plist
    
    if length(plist) == 1
        p_id = 1;
    else
        p_id = id;
        disp(id);
    end
    
    % Reset extrap flag to initial value
    extrap  = extrap_flag;
    
    [~,istop] = min( abs( time - max(t_stop_array) ) );
%     istop = find( time > max(t_stop_array), 1, 'first' );
    
    if temp_final >= temp(istop,p_id)
        if ~isempty(find( temp(:,p_id) >= temp_final, 1, 'last' ))
            istop = find( temp(:,p_id) >= temp_final, 1, 'last' );
        end
        t_stop   = time(istop);
        extrap   = false;
    elseif temp(istop,p_id) >= temp_nse
        t_stop   = time(istop);
        extrap   = false;
    else
        t_stop   = time(istop);
    end
    
    % This is the time which is written to the profile to begin post-processing
    t_start = t_start_array(p_id);
    
    if istop <= i_peak(p_id)
        peak = 1;
    else
        peak = i_peak(p_id);
    end
    
    % This is the time from which we may extrapolate
    t_peak  = time(peak);

    if extrap || plot_flag
        rho = density(:,p_id);
        [~,iwin_min] = min( abs( time - max( t_stop - time_extrap_max, t_peak ) ) );
        span_tau = floor( sgwin_tau * (istop - iwin_min + 1) / ( time(istop) - time(iwin_min) ) );
        if ~mod(span_tau,2)
            span_tau = span_tau + 1;
        end
        span_tau = max( span_tau, sgorder_tau+1 );
        iwin_min = max( iwin_min-span_tau, 1 );
        iwin_max = min( istop+span_tau, length(time) );
        rho(iwin_min:iwin_max) = smooth( time(iwin_min:iwin_max), density(iwin_min:iwin_max,p_id), span_tau, 'sgolay', sgorder_tau );
        rhodot = gradient( rho, time );
        tau_rho = -rho ./ rhodot;

        adiabatic_raw = 0.34 * temp(:,p_id).^3 ./ density(:,p_id);
        adiabatic = adiabatic_raw;
        [~,iwin_min] = min( abs( time - max( t_stop - time_extrap_max, t_peak ) ) );
        span_adi = floor( sgwin_adi * (istop - iwin_min + 1) / ( time(istop) - time(iwin_min) ) );
        if ~mod(span_adi,2)
            span_adi = span_adi + 1;
        end
        span_adi = max( span_adi, sgorder_adi+1 );
        iwin_min = max( iwin_min-span_adi, 1 );
        iwin_max = min( istop+span_adi, length(time) );
        adiabatic(iwin_min:iwin_max) = smooth( time(iwin_min:iwin_max), adiabatic_raw(iwin_min:iwin_max), span_adi, 'sgolay', sgorder_adi );
% %         rhodot = gradient( density(:,p_id), time );
% %         rhodot = [ 0; diff(density(:,p_id)) ./ diff(time) ];
% %         rhodot = deriv_3pt( density(:,p_id), time );
%         [~,span_tau] = min( abs( time - t_stop + sgwin_tau ) );
%         span_tau = max( istop - span_tau + 1, sgorder_tau + 1 )
%         tic;
%         rho = smooth( time, density(:,p_id), span_tau, 'sgolay', sgorder_tau );
%         toc;
%         rhodot = gradient( rho, time );
% %         rhodot = deriv_5pt( rho, time );
%         
% %         tau_rho = -density(:,p_id) ./ rhodot;
%         tau_rho = -rho ./ rhodot;
%         
% %         temp_smooth = smooth( time, temp(:,p_id), sgwin_adi, 'sgolay', sgorder_adi );
% 
%         [~,span_adi] = min( abs( time - t_stop + sgwin_adi ) );
%         span_adi = max( istop - span_adi + 1, sgorder_adi + 1 )
%         
%         adiabatic_raw = 0.34 * temp(:,p_id).^3 ./ density(:,p_id);
% %         adiabatic_raw = temp_smooth .* rho.^(-1/3);
% 
%         tic;
%         adiabatic   = smooth( time, adiabatic_raw, span_adi, 'sgolay', sgorder_adi );
%         toc;
% %         adiabatic   = smooth( adiabatic_raw, 'rlowess', sgwin_adi );
% %         adiabatic   = adiabatic_raw;

    end

    
    if plot_flag
%         tau_rho_avg = smooth( tau_rho, 'rlowess' );
%         tau_rho_avg = smooth( tau_rho, sgwin_tau, 'sgolay', sgorder_tau);
        tau_rho_avg = tau_rho;
        
        fig_h = figure;
        axis_h(j,1) = subplot( 'Position', [.1,.30,.8,.60] );
        axis_h(j,2) = subplot( 'Position', [.1,.10,.8,.20] );
        set( axis_h(j,1), 'YLim', [0 temp_nse], ...
            'XLim', [ t_peak - time_bounce, t_stop - time_bounce], ...
            'XTickLabel', [], ...
            'XMinorTick', 'on', ...
            'XGrid', 'on', ...
            'YGrid', 'on', ...
            'Box', 'on', ...
            'NextPlot', 'add', ...
            'FontSize', 18, ...
            'TickLabelInterpreter_I','latex');
        
        set( axis_h(j,2), 'YLim', [-2.99999, 2.99999], ...
            'YTick', [-2,2], ...
            'XLim', [ t_peak - time_bounce, t_stop - time_bounce ], ...
            'XMinorTick', 'on', ...
            'XGrid', 'on', ...
            'YGrid', 'on', ...
            'Box', 'on', ...
            'NextPlot', 'add', ...
            'FontSize', 18, ...
            'TickLabelInterpreter_I','latex');
        
        if time_bounce ~= 0.0
            xlabel( axis_h(j,2), 'Time after bounce [s]' );
        else
            xlabel( axis_h(j,2), 'Elapsed time [s]' );
        end
        ylabel( axis_h(j,1), 'Temperature [GK]' );
        ylabel( axis_h(j,2), 'Expansion Timescale [s]' );
        
        % Plot the non-extrapolated temperature profile
        plot( axis_h(j,1), time(peak:istop) - time_bounce, temp(peak:istop,p_id), ...
            'Color', 'k', ...
            'LineStyle', '-', ...
            'LineWidth', 1.5, ...
            'DisplayName','Temperature');
        
        % Plot the expansion timescale
        plot( axis_h(j,2), time(peak+5:istop)-time_bounce, tau_rho_avg(peak+5:istop), ...
            'Color', 'k', ...
            'LineStyle', '--', ...
            'DisplayName', '${\tau}_{\mathrm{exp}}$' );
        
        axis_h(j,3) = axes( 'Parent', fig_h, ...
            'Position', get(axis_h(j,1),'Position'), ...
            'Color', 'none', ...
            'YLim', [ min(adiabatic(peak:istop)), max(adiabatic(peak:istop)) ], ...
            'YAxisLocation', 'right', ...
            'XTick', [], ...
            'XTickLabel', [], ...
            'Box', 'off', ...
            'NextPlot', 'add', ...
            'FontSize', 18, ...
            'TickLabelInterpreter_I','latex');
        
        % Plot the isentropic curve
        plot( axis_h(j,3), time(peak+5:istop) - time_bounce, adiabatic(peak+5:istop), ...
            'Color', 'k', ...
            'LineStyle', '-.', ...
            'LineWidth', 0.5, ...
            'DisplayName', '$0.34 T_{9}^{3} / \rho_{5}$' );

        num_ticks = length( get( axis_h(j,1), 'YTickLabel' ) );
        ylim = get( axis_h(j,3), 'YLim' );
        ytick = linspace( ylim(1), ylim(2), num_ticks - 2);
        set( axis_h(j,3), 'YTick', ytick );
        set( axis_h(j,3), 'YLim', [ ylim(1) - (ytick(2)-ytick(1)), ylim(2) + (ytick(end)-ytick(end-1)) ] );
        ylabel( axis_h(j,3), '$0.34 T_{9}^{3} / \rho_{5}$' );
        
        axis_h(j,4) = axes( 'Parent', fig_h, ...
            'Position', get(axis_h(j,1),'Position'), ...
            'Color', 'none', ...
            'XLim', get(axis_h(j,1),'XLim'), ...
            'YLim', get(axis_h(j,1),'YLim'), ...
            'YAxisLocation', 'right', ...
            'XTick', [], ...
            'XTickLabel', [], ...
            'YTick', [], ...
            'YTickLabel', [], ...
            'Box', 'off', ...
            'NextPlot', 'add', ...
            'FontSize', 18, ...
            'TickLabelInterpreter_I','latex');
    end
    
    for n = 1:length(t_stop_array)
        [~,istop] = min( abs( time - t_stop_array(n) ) );
%         istop = find( time > t_stop_array(n), 1, 'first' );
        extrap = extrap_flag;
        extrap_window = [];
        
        if temp_final >= temp(istop,p_id)
            if ~isempty(find( temp(:,p_id) >= temp_final, 1, 'last' ))
                istop = find( temp(:,p_id) >= temp_final, 1, 'last' );
            end
            tau(n,j) = 0.0;
            t_stop   = time(istop);
            t_exp    = t_stop;
            t_extrap = t_stop;
            extrap   = false;
        elseif temp_nse <= temp(istop,p_id)
            tau(n,j) = 0.0;
            t_stop   = time(istop);
            t_exp    = t_stop;
            t_extrap = t_stop;
            extrap   = false;
        else
            t_stop   = time(istop);
        end
        
        if extrap  
%% Determine the time window to be used to estimate the expansion timescale

            % Set the maximum size of the window
            time_extrap_max = time_extrap_max0;
            if t_stop <= t_peak + time_extrap_max
                time_extrap_max = t_stop - t_peak;
                warning('%s%f', 'Time range since last peak is less than time_extrap_max for t_stop=', t_stop);
            elseif t_stop <= t_peak + time_extrap_min
                time_extrap_max = 0.0;
                warning('%s%f', 'Time range since last peak is less than time_extrap_min for t_stop=', t_stop);
            end

            % Set the time index for the end of the time window
            iwin_stop = istop;

            % Determine the time index for the maximum window size
            [~,iwin_max] = min( abs( time(peak:istop) - t_stop + time_extrap_max ) );
            iwin_max = iwin_max + peak - 1;

            % Determine the time index for the minimum window size
            [~,iwin_min] = min( abs( time(iwin_max:istop) - t_stop + time_extrap_min ) );
            iwin_min = iwin_min + iwin_max - 1;

            % Determine how far back from the end of the simulation the constant isentropic/adibatic criteria holds
            iwin_adi = find( abs( abs( adiabatic(iwin_max:istop)/adiabatic(istop) ) - 1.0 ) >= adiabatic_tol, 1, 'last' );
            iwin_adi = iwin_adi + iwin_max - 1;

            % Set the beginning of the window from isentropic/adiabatic criteria if between iwin_max and iwin_min
            if isempty(iwin_adi)
                iwin_start = iwin_max;
            elseif iwin_adi <= iwin_min
                iwin_start = iwin_adi;
            else
                iwin_start = istop;
                warning('Isentropic conditions not met');
            end

            % First consider time ranges where constant isentropic/adiabatic criteria is true
            adi_window = iwin_start:iwin_stop;

            % Restrict the window to time ranges where expansion timescale is within limits
            tau_window = adi_window( tau_rho(adi_window) > tau_rho_min & tau_rho(adi_window) < tau_rho_max );

            pmin           = prctile(tau_rho(tau_window),prc_min);
            pmax           = prctile(tau_rho(tau_window),prc_max);
%             extrap_window  = tau_window( tau_rho(tau_window) >= pmin - 1.5*(pmax-pmin) & tau_rho(tau_window) <= pmax + 1.5*(pmax-pmin) );

%             tau_rho_tol = 2*std(tau_rho(tau_window));
%             tau_deviation = abs( mean(tau_rho(tau_window)) - tau_rho(tau_window) );
%             extrap_window   = tau_window( tau_deviation <= tau_rho_tol );
            

            % Find any gaps in the time window
            ijump       = [ 0, find( diff(tau_window) > 1 ),  length(tau_window) ];
            ijump_min   = find( diff(ijump) > min_jumpdiff );

            % Exclude extrema from the time window
            for i = ijump_min
                jstart        = ijump(i)+1;
                jstop         = ijump(i+1);
                jump_window   = tau_window(jstart:jstop);
%                 pmin           = prctile(tau_rho(jump_window),prc_min);
%                 pmax           = prctile(tau_rho(jump_window),prc_max);
                extrap_window  = [ jump_window( tau_rho(jump_window) >= pmin - 1.5*(pmax-pmin) & tau_rho(jump_window) <= pmax + 1.5*(pmax-pmin) ), extrap_window ];;
% %                 tau_rho_tol   = 2*std(tau_rho(jump_window));
% %                 tau_deviation = abs( mean(tau_rho(jump_window)) - tau_rho(jump_window) );
% %                 extrap_window   = [ jump_window( tau_deviation <= tau_rho_tol ), extrap_window ];
            end

            if isempty(extrap_window)
                warning('Timescale over isentropic region not suitable for extrapolation');
                tau(n,j) = 0.0;
                t_exp    = t_stop;
                t_extrap = t_stop;
            else
                % Average instantaneous expansion timescales over window to get tau
                tau(n,j) = mean(tau_rho(extrap_window));
                
                % Find the time corresponding to the final temperature
                t_exp    = t_stop - 3*tau(n,j)*log( temp_final / temp(istop,p_id) );
                
                % Extend extrapolation to fixed time if specified
                if t_exp < time_final
                    t_extrap = time_final;
                else
                    t_extrap = t_exp;
                end
            end
        end
        
        if plot_flag || write_flag
            rtemp_old    = 0.0;
            rdensity_old = 0.0;
            rnu_flx_old  = 0.0;
            maxchange    = zeros(istop,1);
            out_mask     = false(istop,1);

%             [~,nu_tmp,~] = fd_fit(nu_e0(:,:,p_id),nu_e1(:,:,p_id),nu_e2(:,:,p_id),'Tol',1e-5,'MaxIter',10);

            t_step      = -3.0 * tau(n,j) * log( 1.0 - change_min );
            if t_step > 0.0
                k_stop        = 1 + int32( ( t_exp - t_stop ) / t_step );
                time_exp      = linspace( t_stop, t_stop + t_step*double(k_stop), k_stop + 1 )';
                
                exp_factor    = exp( -( time_exp - t_stop ) / (3*tau(n,j)) );
                
                % Switch to a power-law extrapolation to get from t_exp to time_final ( T(t) = T_0 (t/t_0)^(-2/3) )
                if time_exp(end) < time_final
                    ratio     = ( 1.0 - change_min ) ^ ( -3.0 / 2.0 );
                    n_stop    = ceil( log( 1.0 + ( time_final - time_exp(end) ) * ( ratio - 1.0 ) / ( ratio * t_step ) ) / log( ratio ) );
                    ti_step   = ( ratio * t_step * ( 1.0 - ratio .^ ( 1:n_stop ) ) / ( 1.0 - ratio ) )';
                    
                    time_extrap = [ time_exp; time_exp(end) + ti_step ];
                    extrap_factor = [ exp_factor; exp_factor(end) * ( 1.0 + ti_step / time_exp(end) ) .^ ( -2.0 / 3.0 ) ];
                else
                    time_extrap = time_exp;
                    extrap_factor = exp_factor;
                end
                
                % Do not let temperature fall below temp_min
                if temp(istop,p_id) * extrap_factor(end) < temp_min
                    kmin          = find( temp(istop,p_id) * extrap_factor >= temp_min, 1, 'last' );
                    time_extrap   = [ time_extrap(1:kmin); time_extrap(end) ];
                    extrap_factor = [ extrap_factor(1:kmin); extrap_factor(kmin) ];
                end
                
                % Use constant temperature/density to get to t_extend
                if t_extend > time_extrap(end)
                    time_extrap(end+1) = t_extend;
                    extrap_factor(end+1) = extrap_factor(end);
                end
                
%                 if t_extrap < time_final
%                     tt_step = -tau(n,j) * log( ( 1.0 - change_min ) * 0.1 );
%                 end
%                 
%                 if t_extend > time_extrap(end)
%                     time_extrap(end+1) = t_extend;
%                     extrap_factor(end+1)  = extrap_factor(end);
%                 end
            else
                k_stop      = 0;
                t_extrap    = t_stop;
                time_extrap = t_stop;
                extrap_factor  = 1.0;
            end
                         
            % Question: Would it be better to use an average v_rad?
            rad_extrap  = radius(istop,p_id) + v_rad(istop,p_id) * ( time_extrap - t_stop );
            r2_factor   = ( radius(istop,p_id) ./ rad_extrap ).^2;
            
            if output_tdel > 0.0
                
                time_output = (t_start:output_tdel:t_stop)';
                
%                 [~,iout] = histc( time_output, time );
%                 iout(iout == 0) = 1;
                
%                 iout = zeros(length(time_output),1);
%                 for i = 1:length(time_output)
%                     [~,iout(i)] = min(abs(time-time_output(i)));
%                 end
                
%                 out_mask(iout)    = true;
%                 out_mask(1)       = true;
%                 out_mask(peak)    = true;
%                 out_mask(istop)   = true;
                
                temp_output = interp1( time, temp(:,p_id), time_output, 'linear', temp(istop,p_id) );
                rho_output  = interp1( time, density(:,p_id), time_output, 'linear', density(istop,p_id) );
                ye_output   = interp1( time, ye(:,p_id), time_output, 'linear', ye(istop,p_id) );
                
                time_out = [ time_output; ...
                             time_extrap ];
                temp_out = [ temp_output; ...
                             temp(istop,p_id) * extrap_factor ];
                rho_out  = [ rho_output; ...
                             density(istop,p_id) * extrap_factor.^3 ];
                ye_out   = [ ye_output; ...
                             repmat( ye(istop,p_id), length(time_extrap), 1 ) ];
                
                nu1_flx_output = interp1( time, flxtot(:,1,p_id), time_output, 'linear', flxtot(istop,1,p_id) );
                nu2_flx_output = interp1( time, flxtot(:,2,p_id), time_output, 'linear', flxtot(istop,2,p_id) );
                nu3_flx_output = interp1( time, flxtot(:,3,p_id), time_output, 'linear', flxtot(istop,3,p_id) );
                nu4_flx_output = interp1( time, flxtot(:,4,p_id), time_output, 'linear', flxtot(istop,4,p_id));

                nu1_flx_out = [ nu1_flx_output; ...
                                flxtot(istop,1,p_id) .* r2_factor ];
                nu2_flx_out = [ nu2_flx_output; ...
                                flxtot(istop,2,p_id) .* r2_factor ];
                nu3_flx_out = [ nu3_flx_output; ...
                                flxtot(istop,3,p_id) .* r2_factor ];
                nu4_flx_out = [ nu4_flx_output; ...
                                flxtot(istop,4,p_id) .* r2_factor ];
                
                % Set fluxes to zero after some time corresponding to the end of neutrino emission
                nu_mask                = ( time_out > nu_time_stop );
                nu1_flx_out( nu_mask ) = 0.0;
                nu2_flx_out( nu_mask ) = 0.0;
                nu3_flx_out( nu_mask ) = 0.0;
                nu4_flx_out( nu_mask ) = 0.0;
                
                nu1_temp_output = interp1( time, nu_temp(:,1,p_id), time_output, 'linear', nu_temp(istop,1,p_id) );
                nu2_temp_output = interp1( time, nu_temp(:,2,p_id), time_output, 'linear', nu_temp(istop,2,p_id) );
                nu3_temp_output = interp1( time, nu_temp(:,3,p_id), time_output, 'linear', nu_temp(istop,3,p_id) );
                nu4_temp_output = interp1( time, nu_temp(:,4,p_id), time_output, 'linear', nu_temp(istop,4,p_id) );
                            
                nu1_temp_out = [ nu1_temp_output; ...
                                 repmat( nu_temp(istop,1,p_id), length(time_extrap), 1 ) ];
                nu2_temp_out = [ nu2_temp_output; ...
                                 repmat( nu_temp(istop,2,p_id), length(time_extrap), 1 ) ];
                nu3_temp_out = [ nu3_temp_output; ...
                                 repmat( nu_temp(istop,3,p_id), length(time_extrap), 1 ) ];
                nu4_temp_out = [ nu4_temp_output; ...
                                 repmat( nu_temp(istop,4,p_id), length(time_extrap), 1 ) ];

            else
                for i = 1:istop
                    change_t = abs( 1.0 - temp(i,p_id) * rtemp_old );
                    change_d = abs( 1.0 - density(i,p_id) * rdensity_old ) * 0.1;
                    change_f = abs( 1.0 - flxtot(i,1,p_id) * rnu_flx_old ) * 0.1;
                    change = [ change_t, change_d, change_f ];
                    maxchange(i) = max(change);
                    if maxchange(i) >= change_min || i == 1 || i == istop
                        rtemp_old    = 1/temp(i,p_id);
                        rdensity_old = 1/density(i,p_id);
                        rnu_flx_old  = 1/flxtot(i,1,p_id);
                    end
                end
            
                out_mask        = ( maxchange >= change_min );
                out_mask(1)     = true;
                out_mask(peak)  = true;
                out_mask(istop) = true;

                time_out = [ time( out_mask ); ...
                             time_extrap ];
                temp_out = [ temp( out_mask, p_id ); ...
                             temp(istop,p_id) * extrap_factor ];
                rho_out  = [ density( out_mask, p_id ); ...
                             density(istop,p_id) * extrap_factor.^3 ];
                ye_out   = [ ye( out_mask, p_id ); ...
                             repmat( ye(istop,p_id), length(time_extrap), 1 ) ];

                nu1_flx_out = [ flxtot( out_mask, 1, p_id ); ...
                                flxtot(istop,1,p_id) .* r2_factor ];
                nu2_flx_out = [ flxtot( out_mask, 2, p_id ); ...
                                flxtot(istop,2,p_id) .* r2_factor ];
                nu3_flx_out = [ flxtot( out_mask, 3, p_id ); ...
                                flxtot(istop,3,p_id) .* r2_factor ];
                nu4_flx_out = [ flxtot( out_mask, 4, p_id ); ...
                                flxtot(istop,4,p_id) .* r2_factor ];
                
                % Set fluxes to zero after some time corresponding to the end of neutrino emission
                nu_mask                = ( time_out > nu_time_stop );
                nu1_flx_out( nu_mask ) = 0.0;
                nu2_flx_out( nu_mask ) = 0.0;
                nu3_flx_out( nu_mask ) = 0.0;
                nu4_flx_out( nu_mask ) = 0.0;

                nu1_temp_out = [ nu_temp( out_mask, 1, p_id ); ...
                                 repmat( nu_temp(istop,1,p_id), length(time_extrap), 1 ) ];
                nu2_temp_out = [ nu_temp( out_mask, 2, p_id ); ...
                                 repmat( nu_temp(istop,2,p_id), length(time_extrap), 1 ) ];
                nu3_temp_out = [ nu_temp( out_mask, 3, p_id ); ...
                                 repmat( nu_temp(istop,3,p_id), length(time_extrap), 1 ) ];
                nu4_temp_out = [ nu_temp( out_mask, 4, p_id ); ...
                                 repmat( nu_temp(istop,4,p_id), length(time_extrap), 1 ) ];
                         
            end

            th_out = [ time_out'; ...
                temp_out'; ...
                rho_out'; ...
                ye_out'; ...
                nu1_flx_out'; ...
                nu2_flx_out'; ...
                nu3_flx_out'; ...
                nu4_flx_out'; ...
                nu1_temp_out'; ...
                nu2_temp_out'; ...
                nu3_temp_out'; ...
                nu4_temp_out' ];
            
            [~,ia,~] = unique( th_out(1,:) );
            th_out   = th_out(:,ia);
        end
        
        if write_flag
            
            t_write  = th_out(1,end);
            
            th_fname = sprintf('%s%05.0f%s%04.0f%s', th_fname_base, id, '_', t_write*1000, 'ms');
            
            f_id = fopen( th_fname, 'w' );
            fprintf( f_id, '%s%05.0f\n', 'Tracer Particle ', id );
            fprintf( f_id, '% 15.7E%s\n', t_start, ' Start Time' );
            fprintf( f_id, '% 15.7E%s\n', t_write, ' Stop Time'  );
            fprintf( f_id, '% 15.7E%s\n', 1e-6*(t_write-t_start), ' Initial Timestep' );
            i_outstart = find( th_out(1,:) >= t_start, 1, 'first' );
            i_outstop  = find( th_out(1,:) >= t_write, 1, 'first' );
            fprintf( f_id, th_format, th_out(:,i_outstart:i_outstop) );
            fclose(f_id);
            
        end
        
        if plot_flag && ~isempty(extrap_window)
            
            % Highlight the time range used for extrapolation on plot
%             plot( axis_h(j,4), time(extrap_window) - time_bounce, zeros(size(extrap_window)), ...
%                 'Color', color_array(n,:), ...
%                 'Marker', 'o', ...
%                 'MarkerSize', 4, ...
%                 'MarkerEdgeColor', 'none', ...
%                 'MarkerFaceColor', color_array(n,:), ...
%                 'LineStyle', 'none', ...
%                 'DisplayName', 'Range for extrap' );
            plot( axis_h(j,2), time(extrap_window) - time_bounce, tau_rho_avg(extrap_window), ...
                'Color', color_array(n,:), ...
                'Marker', '.', ...
                'LineStyle', 'none', ...
                'DisplayName', sprintf( '%s%1.4f%s', '${{\tau}^{*}}_{\mathrm{exp}}(t_{f} =$ ', t_stop - time_bounce, ' s)' ) );
            
            % Plot the extrapolated temperature profile
            plot( axis_h(j,1), time_out(time_out>time(istop)) - time_bounce, temp_out(time_out>time(istop)), ...
                'Color', color_array(n,:), ...
                'LineWidth', 1.5, ...
                'DisplayName', sprintf( '%s%1.4f%s', '$T(t; <{{\tau}^{*}}_{\mathrm{exp}}> =$ ', tau(n,j), ' s)' ) );
            
            if t_stop > time(min(extrap_window))
                plot( axis_h(j,3), [ time(min(extrap_window)), t_stop ] - time_bounce, repmat( adiabatic(istop), 1, 2 ), ...
                    'Marker', 'o', ...
                    'MarkerSize', 6, ...
                    'MarkerFaceColor', 'none', ...
                    'MarkerEdgeColor', color_array(n,:), ...
                    'LineWidth', 0.5, ...
                    'LineStyle', '--', ...
                    'Color', color_array(n,:), ...
                    'DisplayName', adiabatic_string );
            else
                plot( axis_h(j,3), [ time(min(extrap_window)); t_stop ] - time_bounce, adiabatic(istop), ...
                    'Marker', 'o', ...
                    'MarkerSize', 6, ...
                    'MarkerFaceColor', 'none', ...
                    'MarkerEdgeColor', color_array(n,:), ...
                    'LineWidth', 0.5, ...
                    'LineStyle', '--', ...
                    'Color', color_array(n,:), ...
                    'DisplayName', adiabatic_string );
            end
                
        end
    end
    if plot_flag
%         title_string = sprintf( '%s%d%s%1.4f%s%s%1.4f', 'Particle ', p_id, ...
%             '     <{{\tau}^{*}}_{exp}> = ', tau(n,j), ' s', ...
%             '     Y_{e} @ NSE: ', ye_nse(p_id) );
        title_string = sprintf( '%s%d', 'Particle ', id );
        title(title_string);
        
        set( axis_h(j,:), 'XLim', [ t_peak - time_bounce, min( time_out(end), t_stop_array(end) + 0.5 ) ] );
        
        leg_h(1) = legend( axis_h(j,1), 'Location', 'East' );
        leg_h(2) = legend( axis_h(j,2), 'Location', 'East' );
        leg_h(3) = legend( axis_h(j,3), 'Location', 'NorthEast' );
        
        set( leg_h, 'Interpreter', 'latex', ...
                    'Color', 'none', ...
                    'EdgeColor', 'none', ...
                    'FontSize', 18 );
        
        linkaxes( axis_h(j,:), 'x' );
        
        set( axis_h(j,3:4), 'Position', get( axis_h(j,1), 'Position' ) );
%         set( axis_h(j,4), 'Position', get( axis_h(j,2), 'Position' ) );
        
        if print_flag
            plot_fname  = sprintf( '%s%s%d%s', plot_folder, 'p', id, '-T9_extrapolation' );
            print( '-dpdf', plot_fname );
        end
    end
    
    j = j + 1;
end

%% Finalize
varargout{1} = tau;
if nargout > 1
    if write_flag
        varargout{2} = th_out(:,i_outstart:i_outstop);
    else
        varargout{2} = th_out;
    end
end

end