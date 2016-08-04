%% find_expansion_timescale
% find_expansion_timescale generates expansion timescales for thermodynamic extrapolations
function [ tau ] = find_expansion_timescale( time, temp, density, plist, varargin )
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
%   plist           Integer      N          1D array containing particle IDs to be used
%     >0                                        (eg. 1:4000 or 1:40:4000 or [70,543,800,...]
%
%% Output Arguments
%   Variable        Type        Dimension   Description
%   ---------       -------     ---------   -------------------------------------------------
%   tau             Float        N          Expansion timescale for each particle in 'plist' [s]
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
%   ExtrapTemp      temp_final      Float        1          Final temperature in extrapolation [GK]
%     >0                                                        ( Default: 0.5 )
%
%   ExtrapTime      time_final      Float        1          Final time in extrapolation [s]
%     >0                                                        ( Default: time to reach temp_final )
%
%   TempNSE         temp_nse        Float        1          Temperature from which to perform extrapolation [GK]
%     >0                                                        ( Default: 8.0 )
%
%   TimeStart       t_start_array   Float        N or 1     Initial time for post-processing profile. [s]
%     >=0                                                       If N=1, the same value is used for all particles
%                                                               ( Default: last time when temp < temp_nse )
%
%   TimeStop        t_stop_array    Float        N          Times from which to perform extrapolation
%                                                               ( Default: time(end) )
%
%% Initialization

% Create an instance of the inputParser class.
p = inputParser;
p.FunctionName = 'find_expansion_timescale';
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
p.addRequired('plist', ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'vector', ...
                                         'integer', ...
                                         '>', 0}));
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
p.addOptional('TempNSE', 8.0, ...
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

% Parse, validate, and assign required input arguments.
p.parse( time, temp, density, plist, varargin{:} );
time            = p.Results.time;
temp            = p.Results.temp;
density         = p.Results.density;
plist           = p.Results.plist;
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
temp_nse        = p.Results.TempNSE;
temp_final      = p.Results.ExtrapTemp;
time_final      = p.Results.ExtrapTime;
t_start_array   = p.Results.TimeStart;
t_stop_array    = p.Results.TimeStop;

prc_min = 25; prc_max = 75;

time_extrap_max0 = time_extrap_max;

if t_stop_array == 0.0
    t_stop_array = time(end);
end

if time_final > 0.0
    t_stop_array(t_stop_array > time_final) = time_final;
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
        [~,i_nse]    = min( abs( time - t_start_array(p_id) ) );
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

j = 1;
for id = plist
    
    if length(plist) == 1
        p_id = 1;
    else
        p_id = id;
        disp(id);
    end
    
    istop_array(1:length(t_stop_array)) = length(time);
    for n = 1:length(t_stop_array)
        [~,istop_array(n)] = min( abs( time - t_stop_array(n) ) );
    end
    
    if all( temp(istop_array,p_id) >= temp_nse ) || all( temp_final >= temp(istop_array,p_id) )
        tau(:,j) = 0.0;
        continue
    end
    
    % This is the time from which we may extrapolate
%     [~,istop] = min( abs( time - max(t_stop_array) ) );
    istop = max(istop_array);
%     istop = find( time > max(t_stop_array), 1, 'first' );
    
    if temp_final >= temp(istop,p_id)
        if ~isempty(find( temp(:,p_id) >= temp_final, 1, 'last' ))
            istop = find( temp(:,p_id) >= temp_final, 1, 'last' );
        end
        t_stop   = time(istop);
    elseif temp(istop,p_id) >= temp_nse
        t_stop   = time(istop);
    else
        t_stop   = time(istop);
    end
%     
%     % This is the time which is written to the profile to begin post-processing
%     t_start = t_start_array(p_id);
%     
    if istop <= i_peak(p_id)
        peak = 1;
    else
        peak = i_peak(p_id);
    end
    
    % This is the time from which we may extrapolate
    t_peak  = time(peak);
    
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
    
    for n = 1:length(t_stop_array)
%% Determine the time window to be used to estimate the expansion timescale
        
        % Determine the time from which to perform the extrapolation
%         [~,istop] = min( abs( time - t_stop_array(n) ) );
        istop = istop_array(n);
        extrap_window = [];
        
        if temp_final >= temp(istop,p_id)
            if ~isempty(find( temp(:,p_id) >= temp_final, 1, 'last' ))
                istop = find( temp(:,p_id) >= temp_final, 1, 'last' );
            end
            tau(n,j) = 0.0;
            t_stop   = time(istop);
        elseif temp_nse <= temp(istop,p_id)
            tau(n,j) = 0.0;
            t_stop   = time(istop);
        else
            t_stop   = time(istop);
        end
        
        % Set the maximum size of the window
        time_extrap_max = time_extrap_max0;
        if t_stop <= t_peak + time_extrap_max
            time_extrap_max = t_stop - t_peak;
%             warning('%s%f', 'Time range since last peak is less than time_extrap_max for t_stop=', t_stop);
        elseif t_stop <= t_peak + time_extrap_min
            time_extrap_max = 0.0;
%             warning('%s%f', 'Time range since last peak is less than time_extrap_min for t_stop=', t_stop);
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
%             warning('Isentropic conditions not met');
        end
        
        % First consider time ranges where constant isentropic/adiabatic criteria is true
        adi_window = iwin_start:iwin_stop;
        
        % Restrict the window to time ranges where expansion timescale is within limits
        tau_window = adi_window( tau_rho(adi_window) > tau_rho_min & tau_rho(adi_window) < tau_rho_max );
        
        pmin           = prctile(tau_rho(tau_window),prc_min);
        pmax           = prctile(tau_rho(tau_window),prc_max);
%         extrap_window  = tau_window( (tau_rho(tau_window) >= pmin - 1.5*(pmax-pmin)) & (tau_rho(tau_window) <= pmax + 1.5*(pmax-pmin)) );
        
%         tau_rho_tol = 2*std(tau_rho(tau_window));
%         tau_deviation = abs( mean(tau_rho(tau_window)) - tau_rho(tau_window) );
%         extrap_window   = tau_window( tau_deviation <= tau_rho_tol );
        
        % Find any gaps in the time window
        ijump       = [ 0, find( diff(tau_window) > 1 ),  length(tau_window) ];
        ijump_min   = find( diff(ijump) > min_jumpdiff );
        
        % Exclude extrema from the time window
        for i = ijump_min
            jstart         = ijump(i)+1;
            jstop          = ijump(i+1);
            jump_window    = tau_window(jstart:jstop);
%             pmin           = prctile(tau_rho(jump_window),prc_min);
%             pmax           = prctile(tau_rho(jump_window),prc_max);
            extrap_window  = [ jump_window( tau_rho(jump_window) >= pmin - 1.5*(pmax-pmin) & tau_rho(jump_window) <= pmax + 1.5*(pmax-pmin) ), extrap_window ];
% %             tau_rho_tol    = 2*std(tau_rho(jump_window));
% %             tau_deviation  = abs( mean(tau_rho(jump_window)) - tau_rho(jump_window) );
%             extrap_window  = [ jump_window( tau_deviation <= tau_rho_tol ), extrap_window ];
        end
        
        if isempty(extrap_window)
%             warning('Timescale over isentropic region not suitable for extrapolation');
            tau(n,j) = 0.0;
        else
            tau(n,j) = mean(tau_rho(extrap_window));
        end
        
    end
    
    j = j + 1;
    
end

