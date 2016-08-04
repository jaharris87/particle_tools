function [ varargout ] = plot_fate_bins( time, p_negvr, p_posvr, p_bound, p_pns, fate_bin, varargin )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Create an instance of the inputParser class.
p = inputParser;

% Define required inputs
p.addRequired('time', ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'vector', ...
                                         'real', ...
                                         '>=', 0}));
p.addRequired('p_negvr', ...
              @(x)validateattributes(x, {'cell'}, ...
                                        {'vector'}));
p.addRequired('p_posvr', ...
              @(x)validateattributes(x, {'cell'}, ...
                                        {'vector'}));
p.addRequired('p_bound', ...
              @(x)validateattributes(x, {'cell'}, ...
                                        {'vector'}));
p.addRequired('p_pns', ...
              @(x)validateattributes(x, {'cell'}, ...
                                        {'vector'}));
p.addRequired('fate_bin', ...
              @(x)validateattributes(x, {'char'}, ...
                                        {'row'}));

% Parse, validate, and assign required input arguments.
p.parse( time, p_negvr, p_posvr, p_bound, p_pns, fate_bin );
time     = p.Results.time;
p_negvr  = p.Results.p_negvr;
p_posvr  = p.Results.p_posvr;
p_bound  = p.Results.p_bound;
p_pns    = p.Results.p_pns;
fate_bin = p.Results.fate_bin;

switch fate_bin
    case 'bound'
        p_bin = p_bound;
    case 'posvr'
        p_bin = p_posvr;
    case 'negvr'
        p_bin = p_negvr;
    case 'pns'
        p_bin = p_pns;
    otherwise
        error( 'Unexpected fate_bin value. Expected values are "bound", "posvr", "negvr", or "pns"' );
end

ymin = 0;
ymax = max( cellfun( @length, p_bin ) );
ypow10 = floor( log10( ymax ) );
ymax = roundn( ymax + 10^(ypow10-1), ypow10-1 );

                                    
% Define optional inputs
p.addOptional('TimeStart', min(time), ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'real', ...
                                         'nonnegative', ...
                                         '>=', min(time), ...
                                         '<', max(time)}));
p.addOptional('TimeStop', max(time), ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'real', ...
                                         'nonnegative', ...
                                         '>', min(time), ...
                                         '<=', max(time)}));
p.addOptional('XLim', [min(time), max(time)], ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'vector', ...
                                        'real', ...
                                        'numel', 2}));
p.addOptional('YLim', [ymin+1e-5, ymax-1e-5 ], ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'vector', ...
                                        'real', ...
                                        'nonnegative', ...
                                        'numel', 2}));
p.addOptional('TimeWindow', 0.0, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'real'}));
p.addOptional('TimeBounce', 0.0, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'real', ...
                                         'nonnegative'}));
p.addOptional('Normalize', false, ...
              @(x)validateattributes(x, {'logical'}, ...
                                        {'scalar'}));
p.addOptional('PlotTitle', '', ...
              @(x)validateattributes(x, {'char'}, ...
                                        {'row'}));
p.addOptional('TracerMass', 1.0, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'real', ...
                                         'nonnegative'}));
p.addOptional('CreatePlot', true, ...
              @(x)validateattributes(x, {'logical'}, ...
                                        {'scalar'}));

% Parse, validate, and assign optional input arguments
p.parse( time, p_negvr, p_posvr, p_bound, p_pns, fate_bin, varargin{:} );
tstart     = p.Results.TimeStart;
tstop      = p.Results.TimeStop;
time_bounce    = p.Results.TimeBounce;
xlim       = p.Results.XLim;
ylim       = p.Results.YLim;
twindow    = p.Results.TimeWindow;
normalize  = p.Results.Normalize;
plot_title = p.Results.PlotTitle;
tracer_mass= p.Results.TracerMass;
plot_flag  = p.Results.CreatePlot;

if twindow == 0.0
    t_string = '$t = t_{\mathrm{f}}$';
else
    t_string = '$t$';
end

switch fate_bin
    case 'bound'
        bin_string = strcat( '${{\bf{P}}_{\mathrm{bound}}}$(', t_string, ')' );
    case 'posvr'
        bin_string = strcat( '${{\bf{P}}_{\mathrm{unb}}^{+}}$(', t_string, ')' );
    case 'negvr'
        bin_string = strcat( '${{\bf{P}}_{\mathrm{unb}}^{-}}$(', t_string, ')' );
    case 'pns'
        bin_string = strcat( '${{\bf{P}}_{\mathrm{PNS}}^{-}}$(', t_string, ')' );
    otherwise
        error( 'Unexpected fate_bin value. Expected values are "bound", "posvr", "negvr", or "pns"' );
end

load_constants;

% Calculate indices for plotting
[~,istart] = min( abs( time - tstart ) );
[~,istop]  = min( abs( time - tstop ) );

xplot = time(istart:istop) - time_bounce;
yplot = zeros(length(xplot),5);

if twindow == 0.0
    p_bin(1:istop) = p_bin(istop);
end

if normalize
    norm_factor = 100.0 ./ cellfun( @length, p_bin );
    ylim = [1e-5, 100-1e-5];
else
    norm_factor = ones( 1, length(time) );
end

[~,bounce_index] = min( abs( time - time_bounce ) );
if twindow >= 0.0
    window = 1;
    ii_start = bounce_index+1;
    ii_stop = istop;
    ii_inc  = 1;
else
    window = -1;
    ii_start = istop-1;
    ii_stop  = bounce_index;
    ii_inc   = -1;
end

for i = ii_start:ii_inc:ii_stop
    if twindow > 0.0
        if abs( time(i-window) - time(i) - twindow ) / twindow > 0.001
            [~,itime] = min( abs( time - time(i) + twindow ) );
            window    = i - itime + 1;
        end
    elseif twindow < 0.0
        if abs( time(i-window) - time(i) - twindow ) / twindow < 0.001
            [~,itime] = min( abs( time - time(i) + twindow ) );
            window    = (i - itime + 1);
        end
    end
    j = i - istart + 1;
    yplot(j,1) = length( intersect( p_bin{i}, p_posvr{i-window} ) );
    yplot(j,2) = length( intersect( p_bin{i}, p_negvr{i-window} ) );
    yplot(j,3) = length( intersect( p_bin{i}, p_pns{i-window} ) );
    yplot(j,4) = length( intersect( p_bin{i}, p_bound{i-window} ) );
    yplot(j,5) = length( p_bin{i} );
    yplot(j,:) = yplot(j,:) * norm_factor(i);
end
if tracer_mass ~= 1 && ~normalize
    yplot = yplot * tracer_mass / M_solar;
    ylim = ylim * tracer_mass / M_solar;
end

if plot_flag
    fig_h  = figure( 'PaperOrientation', 'landscape', ...
                     'PaperPosition', [-0.1 0.2 8.3 3.5], ...
                     'PaperSize', [7.65, 3.6], ...
                     'PaperUnits', 'inches' );

    % Set up the axis for the plot
    axis_h(1) = axes( 'Parent', fig_h, ...
                      'XLim', xlim, ...
                      'XMinorTick', 'on', ...
                      'YLim', ylim, ...
                      'YMinorTick', 'on', ...
                      'Box', 'on', ...
                      'NextPlot', 'replacechildren', ...
                      'FontSize', 18, ...
                      'TickLabelInterpreter_I', 'latex', ...
                      'Position', [0.1, 0.15, 0.725, 0.8] );

    if ~isempty(plot_title)
    %     if twindow == 0.0
    %         plot_title = sprintf( '%s%s%s%d%s', 'Fates of ', bin_string, ' Particles at {\tau} = t' );
    %     else
    %         plot_title = sprintf( '%s%s%s%d%s', 'Fates of ', bin_string, ' Particles at {\tau} = t - ', round(twindow*1.0E+03), ' ms' );
    %     end
        title( axis_h(1), plot_title, 'FontSize', 18 );
    end

    % grid( axis_h(1), 'on' );

    if tracer_mass ~= 1
        if normalize
            ylabel( axis_h(1), '\% of Mass', 'FontSize', 18 );
        else
            ylabel( axis_h(1), 'Mass [$\times 10^{-2}$ M$_\odot$]', 'FontSize', 18 );
        end
    else
        if normalize
            ylabel( axis_h(1), '\% of Particles', 'FontSize', 18 );
        else
            ylabel( axis_h(1), '\# of Particles', 'FontSize', 18 );
        end
    end

    if( time_bounce == 0 )
        xlabel( axis_h(1), 'Elapsed time [ms]', 'FontSize', 18 );
    else
        xlabel( axis_h(1), 'Time after bounce [s]', 'FontSize', 18 );
    end
    % ColorSet = varycolor( length( yplot(1,:) ) );
    ColorSet(1,:) = [255 0 0];
    ColorSet(2,:) = [0 0 255];
    ColorSet(3,:) = [237 137 32];
    ColorSet(4,:) = [0 127 0];
    ColorSet(5,:) = [0 0 0];
    ColorSet = ColorSet / 255;

    set( axis_h(1), 'ColorOrder', ColorSet );

    line_h = plot( axis_h(1), xplot, yplot(:,1:5), ...
                              'LineWidth', 2.0, ...
                              'LineStyle', '-' );

    % Add a legend
    if twindow > 0.0
        tau_string = sprintf( '%s%d%s','$t - ', abs(round(twindow*1.0E+03)), '$ ms' );
    elseif twindow < 0.0
        tau_string = sprintf( '%s%d%s','$t + ', abs(round(twindow*1.0E+03)), '$ ms' );
    else
        tau_string = '$t$';
    end
    posvr_label = strcat( '${{\bf{P}}_{\mathrm{unb}}^{+}}$(', tau_string, ')' );
    negvr_label = strcat( '${{\bf{P}}_{\mathrm{unb}}^{-}}$(', tau_string, ')' );
    pns_label   = strcat( '${{\bf{P}}_{\mathrm{PNS}}}$(', tau_string, ')' );
    bound_label = strcat( '${{\bf{P}}_{\mathrm{bound}}}$(', tau_string, ')' );
    leg_h = legend( line_h, {posvr_label; ...
                             negvr_label; ...
                             pns_label; ...
                             bound_label; ...
                             bin_string} );
    set( leg_h, 'Location', 'NorthWest', ...
                'Interpreter', 'latex', ...
                'Color', 'none', ...
                'EdgeColor', 'none', ...
                'FontSize', 18 );
end

% Output the axis and line object handles
varargout{1} = yplot;
varargout{2} = line_h;
varargout{3} = leg_h;
                         
end

