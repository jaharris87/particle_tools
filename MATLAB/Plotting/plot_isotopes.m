%% plot_isotopes
% Plot isotope quantities (e.g. mass or abundance) for each element (i.e. Z) as a function of mass number (i.e. A)
function [ varargout ] = plot_isotopes( aa, zz, mx, varargin )
%% Input Arguments
%   Variable        Type        Dimension   Units   Description
%   --------        -------     ---------   -----   ------------------------------------------------------------
%   aa              Integer      M                  Mass number (A) for each isotope
%     >0
% 
%   zz              Integer      M                  Atomic/proton number (Z) for each isotope
%     >=0
%
%   mx              Float       (M, N)              Each column to be separately plotted for 1 <= N <= 4
%     >=0.0
%
%% Output Arguments
%   Variable        Type        Dimension   Units   Description
%   ---------       -------     ---------   -----   ------------------------------------------------------------
%   zz_plot         Cell         K                  List of unique Z values
%
%   aa_plot         Cell         K                  Each cell element i contains values of A for each Z
%
%   mx_plot         Cell        (K, N)              Each cell element (i,j) contains values of mx for each Z
%
%   saxis_h         Float        L                  Axis handles for each sub-axis
%
%   text_h          Cell         L                  Text-box object handles for element names for each Z
%
%   mx_h            Cell        (N, L)              Line object handles for each Z
%
%% Optional Input Arguments
%   Flag            Variable        Type        Dimension   Units   Description                                                     Default
%   ---------       --------        -------     ---------   -----   ------------------------------------------------------------    ----------------------------
%   ALim            aa_lim          Integer      2                  Min and max values of A to plot                                 [0, max(aa)];
%    >0
%
%   MinMx           mx_min          Float        1                  Values in mx smaller than mx_min will be marked and plotted     0
%    >=0.0                                                              separately
%
%   Delta           delta           Float        1          norm    Vertical spacing between subplots                               0.08
%    >=0.0
%
%   YLim            ylim            Float        2                  Axis limits for y-axis                                          [min(mx_plot),max(mx_plot)];
%    >0
%
%   YScale          yscale          String       1                  Scale to use for y-axis                                         'log'
%    'linear'
%    'log'
%
%% Initialization

load_constants;

% Create an instance of the inputParser class.
p = inputParser;
p.FunctionName = 'plot_isotopes';
p.KeepUnmatched = true;

% Define required inputs
p.addRequired('aa', ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'vector', ...
                                         'integer', ...
                                         '>', 0}));
p.addRequired('zz', ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'vector', ...
                                         'integer', ...
                                         '>=', 0}));
p.addRequired('mx', ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'2d', ...
                                         'real'}));
% Define optional inputs
p.addOptional('ALim', [], ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'vector', ...
                                         'integer', ...
                                         'numel', 2, ...
                                         '>=', 0}));
p.addOptional('MinMx', 0.0, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'real', ...
                                         '>=', -Inf}));
p.addOptional('Delta', 0.08, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'real', ...
                                         '>=', 0.0}));
p.addOptional('XOverlap', 0, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'integer', ...
                                         '>=', 0}));
p.addOptional('XWidth', 0, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'integer', ...
                                         '>=', 0}));
p.addOptional('XTickSpacing', 2, ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'scalar', ...
                                         'integer', ...
                                         '>', 0}));
p.addOptional('YLim', [], ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'vector', ...
                                         'real', ...
                                         'numel', 2}));
p.addOptional('YScale', 'log', ...
              @(x)validateattributes(x, {'char'}, ...
                                        {'row'}));
p.addOptional('Solar', false, ...
              @(x)validateattributes(x, {'logical'}, ...
                                        {'scalar'}));
p.addOptional('LogEps', false, ...
              @(x)validateattributes(x, {'logical'}, ...
                                        {'scalar'}));

% Parse, validate, and assign input arguments.
p.parse( aa, zz, mx, varargin{:} );
aa = p.Results.aa;
zz = p.Results.zz;
mx = p.Results.mx;

aa_lim   = p.Results.ALim;
mx_min   = p.Results.MinMx;
delta    = p.Results.Delta;
ylim     = p.Results.YLim;
yscale   = p.Results.YScale;
solar    = p.Results.Solar;
logeps   = p.Results.LogEps;

xoverlap      = p.Results.XOverlap;
xwidth        = p.Results.XWidth;
xtick_spacing = p.Results.XTickSpacing;

if isempty(aa_lim)
    aa_lim(1) = 0;
    aa_lim(2) = max(aa);
end

% Get the dimensions from input
nmx = size(mx,2);

% Identify unique elements
zz_plot = unique(zz)';
nzz     = length(zz_plot);

% Element label text box vertical alignment
valign = 'middle';

% Plot window position
pos = [ 0.1, 0.1, 0.85, 0.85 ];

% ylim = [ .10001E-12, 9.9999 ];

% Combined domain of horizontal axes
xmin        = aa_lim(1);
xmax        = aa_lim(2);

% Width of each horizontal axis
if xwidth == 0
    xwidth = xmax - xmin;
end

% Number of subplots necessary for horizontal axis parameters
nsubplot = ceil( (xmax-xmin-xoverlap) / (xwidth-xoverlap) );

% Preallocate
aa_plot =  cell( nzz, 1 );
mx_plot =  cell( nzz, nmx );
mx_h    =  cell( nmx, nsubplot );
text_h  =  cell(   1, nsubplot );
saxis_h = zeros(   1, nsubplot );

plot_flag = cell( nzz, 1 );

% Group isotopes of each unique element into cells
for i = 1:nzz
    iz = zz_plot(i);
    
    aa_plot{i} = aa(zz==iz);
    for j = 1:nmx
        mx_plot{i,j} = mx(zz==iz,j);
    end
end

% Vertical axis limits
if isempty(ylim)
    ylim(1) = min( cellfun(@min, mx_plot) );
    ylim(2) = max( cellfun(@max, mx_plot) );
end

% Cell array of element names
element = cellstr( build_element_symbol() );

% Axis labels
xtick         = xmin:xtick_spacing:xmax;
xtick_label   = strsplit(num2str(xtick),' ');
if xtick_spacing <= 2
    xtick_label(2:2:end) = {''};
    xminortick = 'off';
else
    xminortick = 'on';
end

if strcmp(yscale,'log')
    ytick = logspace(floor(log10(ylim(1))),ceil(log10(ylim(2))),ceil(log10(ylim(2)))-floor(log10(ylim(1)))+1);
    ytick_label = strsplit( sprintf( '%g ', ytick ) );
    ytick_label = ytick_label(1:end-1);
end

% Horizontal axes limts
saxis_xmin = xmin:(xwidth-xoverlap):(xmin+(nsubplot-1)*(xwidth-xoverlap))
saxis_xmax = (xmin+xwidth):(xwidth-xoverlap):(xmin+nsubplot*(xwidth-xoverlap)+xoverlap)
xlim(:,1) = saxis_xmin-0.5;
xlim(:,2) = saxis_xmax+0.5;

% Position of each subaxis
saxis_pos(:,1) = repmat( pos(1), nsubplot, 1 );
saxis_pos(:,2) = pos(2) + pos(4) - (1:nsubplot)*pos(4)/nsubplot + delta;
saxis_pos(:,3) = repmat( pos(3), nsubplot, 1 );
saxis_pos(:,4) = repmat( pos(4)/nsubplot-delta, nsubplot, 1 );

% Create a new figure
fig_h = figure;

if solar
    set( fig_h, 'PaperType', 'A4', 'PaperPosition', [ 0.25 0.25 7.767716 11.192913 ] );
    ylabel_text = 'Production Factor $(X/X_{\odot})$';
    leg_loc = 'NorthOutside';
    i_o16 = find(aa==16 & zz==8);
elseif logeps
    set( fig_h, 'PaperType', 'A4', 'PaperPosition', [ 0.25 0.25 7.767716 11.192913 ] );
    ylabel_text = '$\log_{10}(\varepsilon_{x})$';
    leg_loc = 'NorthOutside';
    i_pp = find(aa==1 & zz==1);
else
    ylabel_text = 'Mass [M$_{\odot}$]';
    leg_loc = 'South';
end

% Set up the subplots
for iplot = 1:nsubplot
    saxis_h(iplot) = subplot( nsubplot, 1, iplot, ...
        'Position', saxis_pos(iplot,:), ...
        'NextPlot', 'add', ...
        'XMinorTick', 'off', ...
        'Box', 'on', ...
        'FontSize', 18, ...
        'TickLabelInterpreter_I', 'latex', ...
        'Yscale', yscale, ...
        'YLim', ylim, ...
        'XLim', xlim(iplot,:), ...
        'XTick', xtick, ...
        'XTickLabel', xtick_label, ...
        'XMinorTick', xminortick );
    if strcmp(yscale,'log')
        set( saxis_h(iplot), 'YTick', ytick );
    end
    if solar
        set( saxis_h(iplot), 'YTickLabel', ytick_label );
    end
end

% Set up axis labels
xlabel( saxis_h(nsubplot), 'Mass Number (A)', ...
    'FontSize', 18, ...
    'Interpreter', 'latex' );
y_h = ylabel( saxis_h(1), ylabel_text, ...
    'FontSize', 18, ...
    'Interpreter', 'latex' );

if nsubplot == 4
    set( y_h, 'Units', 'normalized', 'Position', [ -0.06498379786809286 -1.5736434009186056 0 ] );
end

LineSet = { '-'; '--'; '-.'; ':' };
LineWidth = [ 1.5; 1.5; 1.5; 1.5 ];
nlines  = 4;

MarkerSet  = { '.'; 'o'; 'd'; 's' };
MarkerSize = [ 24; 6; 6; 6 ];
nmarkers   = 4;

%ColorSet = get( saxis_h(1), 'ColorOrder' );
ColorSet(1,:)  = [hex2dec('44'),hex2dec('77'),hex2dec('AA')];
ColorSet(2,:)  = [hex2dec('CC'),hex2dec('66'),hex2dec('77')];
ColorSet(3,:)  = [hex2dec('11'),hex2dec('77'),hex2dec('33')];
ColorSet(4,:)  = [hex2dec('DD'),hex2dec('CC'),hex2dec('77')];
ColorSet(5,:)  = [hex2dec('88'),hex2dec('CC'),hex2dec('EE')];
ColorSet(6,:)  = [hex2dec('AA'),hex2dec('44'),hex2dec('99')];
ColorSet(7,:)  = [hex2dec('44'),hex2dec('AA'),hex2dec('99')];
ColorSet(8,:)  = [hex2dec('99'),hex2dec('99'),hex2dec('33')];
ColorSet(9,:)  = [hex2dec('88'),hex2dec('22'),hex2dec('55')];
ColorSet = ColorSet ./ 255;
ncolors  = size( ColorSet, 1 );

% ytext = text(-4.464011669953662,1.9013838345494575E-6,'Mass [M$_{\odot}$]',...
%     'Parent',saxis_h(1),...
%     'Interpreter','latex',...
%     'fontsize',18,...
%     'rotation',90,...
%     'HorizontalAlignment','center',...
%     'verticalalignment','bottom');

legaxis_h = axes( 'Position', saxis_pos(1,:), ...
    'NextPlot', 'add', ...
    'box', 'on', ...
    'XMinorTick', 'off', ...
    'FontSize', 18, ...
    'TickLabelInterpreter_I', 'latex', ...
    'Yscale', yscale, ...
    'YLim', ylim, ...
    'YTick', [], ...
    'YTickLabel', [], ...
    'XLim', xlim(1, :), ...
    'XTick', [], ...
    'XTickLabel', [], ...
    'Color', 'none' );

if nmx == 2
    if solar || logeps
        PlotName = { 'CHIMERA'; ...
                     'WH07' };
    else
        PlotName = { '$\tau^{*}_{\max}$'; ...
                     '$\tau^{*}_{\min}$' };
    end
elseif nmx == 3
    PlotName = { '${{\bf{P}}_{\mathrm{unb}}}$'; ...
                 '${{\bf{P}}_{\mathrm{unb}}^{-}}$'; ...
                 '$\lbrack{{{\bf{P}}_{\mathrm{unb}}^{-}}}\rbrack$' };
end
         
% Generate the legend on a separate axis
jj = 1;
for j = 1:nmx
    if nmx ~= 2 && nmx ~= 3
        PlotName{j,:} = sprintf( '%s%d', 'Plot', j );
    end
    plot( legaxis_h, [-1,-1], [1,1], ...
        'LineStyle', LineSet{jj,:}, ...
        'LineWidth', LineWidth(jj), ...
        'Marker', MarkerSet{jj,:}, ...
        'MarkerSize', MarkerSize(jj), ...
        'MarkerFaceColor', 'auto', ...
        'Color', 'k', ...
        'DisplayName', PlotName{j,:} );
    
    jj = mod( j, nlines ) + 1;
end
leg_h = legend( legaxis_h, 'Orientation', 'horizontal', 'Location', leg_loc );
set( leg_h, 'Interpreter', 'latex', 'Color', 'none', 'EdgeColor', 'none', 'FontSize', 18 );
set( legaxis_h, 'Position', saxis_pos(1,:) );

% Do the plotting
for i = 1:nzz
    iz = zz_plot(i);
    
    plot_flag{i} = false( 1, nsubplot );
    for iplot = 1:nsubplot
        if any( aa_plot{i} <= xlim(iplot,2) ) && any( aa_plot{i} >= xlim(iplot,1) )
            plot_flag{i}(iplot) = true;
        end
    end
end

tmp_flag = true;
k = 1;

xpos = zeros(size(zz_plot));
ypos = NaN(size(zz_plot));
ipos = ones(size(zz_plot));

for i = 1:nzz
    iz = zz_plot(i);
    
    if any( plot_flag{i} )
        
        tmp_ypos = NaN( 1, nmx );
        tmp_ipos = ones( 1, nmx );
        if tmp_flag
%             [ ypos(i), ipos(i) ] = max( max( [mx_plot{i,:}] ) );
            for j = 1:nmx
                mask = mx_plot{i,j} >= mx_min;
                if any( mask )
                    tmp_ypos(j) = max( mx_plot{i,j}(mask) );
                    [ ~, tmp_ipos(j) ] = min( abs( tmp_ypos(j) - mx_plot{i,j} ) );
                else
                    [ tmp_ypos(j), tmp_ipos(j) ] = max( mx_plot{i,j} );
                end
            end
            [ ypos(i), jpos ] = max( tmp_ypos );
            ipos(i) = tmp_ipos(jpos);
        else
%             [ ypos(i), ipos(i) ] = min( min( [mx_plot{i,:}] ) );
            for j = 1:nmx
                mask = mx_plot{i,j} >= mx_min;
                if any( mask )
                    tmp_ypos(j) = min( mx_plot{i,j}(mask) );
                    [ ~, tmp_ipos(j) ] = min( abs( tmp_ypos(j) - mx_plot{i,j} ) );
                else
                    [ tmp_ypos(j), tmp_ipos(j) ] = min( mx_plot{i,j} );
                end
            end
            [ ypos(i), jpos ] = min( tmp_ypos );
            ipos(i) = tmp_ipos(jpos);
        end
        tmp_flag = ~tmp_flag;
        xpos(i) = aa_plot{i}(ipos(i));
        
    end
    
end

% tmp_flag = true;
% 
% for i = 1:nzz
%     iz = zz_plot(i);
%     
%     if any( plot_flag{i} )
%         
%         if tmp_flag
%             if sum( xpos(i) == xpos ) > 1 && ypos(i) < max( ypos( xpos == xpos(i) ) )
%                 if ipos(i) > 1
%                     ipos(i) = ipos(i) - 1;
%                 elseif ipos(i) < length( aa_plot{i} )
%                     ipos(i) = ipos(i) + 1;
%                 end
%                 xpos(i) = aa_plot{i}(ipos(i));
%                 ypos(i) = mx_plot{i,1}(ipos(i),1);
%             end
%         else
%             if sum( xpos(i) == xpos ) > 1 && ypos(i) > min( ypos( xpos == xpos(i) ) )
%                 if ipos(i) > 1
%                     ipos(i) = ipos(i) - 1;
%                 elseif ipos(i) < length( aa_plot{i} )
%                     ipos(i) = ipos(i) + 1;
%                 end
%                 xpos(i) = aa_plot{i}(ipos(i));
%                 ypos(i) = mx_plot{i,1}(ipos(i));
%             end
%         end
%         tmp_flag = ~tmp_flag;
% 
%     end
%     
% end

if solar
    for iplot = 1:nsubplot
        
        o16_h(1) = plot( saxis_h(iplot), xlim(iplot,:), [1 1]*mx(i_o16,1), ...
            'LineStyle', '--', ...
            'LineWidth', 0.5, ...
            'Marker', 'none', ...
            'Color', 'k' );
        o16_h(2) = plot( saxis_h(iplot), xlim(iplot,:), 2*[1 1]*mx(i_o16,1), ...
            'LineStyle', ':', ...
            'LineWidth', 0.5, ...
            'Marker', 'none', ...
            'Color', 'k' );
        o16_h(3) = plot( saxis_h(iplot), xlim(iplot,:), 0.5*[1 1]*mx(i_o16,1), ...
            'LineStyle', ':', ...
            'LineWidth', 0.5, ...
            'Marker', 'none', ...
            'Color', 'k' );
        
        set(get(get(o16_h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        set(get(get(o16_h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        set(get(get(o16_h(3),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end

for iplot = 1:nsubplot
    group_h(iplot,1) = hggroup( 'Parent', saxis_h(iplot), 'DisplayName', '$X < X_{\min}$' );
    group_h(iplot,2) = hggroup( 'Parent', saxis_h(iplot), 'DisplayName', '$X > X_{\min}$' );
%     group_h(iplot,3) = hggroup( 'Parent', saxis_h(iplot) );
    set(get(get(group_h(iplot,1),'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
    set(get(get(group_h(iplot,2),'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
%     set(get(get(group_h(iplot,3),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

tmp_flag = true;

for i = 1:nzz
    iz = zz_plot(i);
    
    if any( plot_flag{i} )
        
%         if tmp_flag
%             
%             if length( mx_plot{i,1} ) > 2;
%                 ipos1 = floor( ( length(mx_plot{i,1}) + 1 ) / 2 );
%                 ypos1 = mx_plot{i,1}(ipos1);
%             else
%                 [ ypos1, ipos1 ] = max( mx_plot{i,1} );
%             end
%             
%             if length( mx_plot{i,2} ) > 2
%                 ipos2 = floor( ( length(mx_plot{i,2}) + 1 ) / 2 );
%                 ypos2 = mx_plot{i,2}(ipos2);
%             else
%                 [ ypos2, ipos2 ] = max( mx_plot{i,2} );
%             end
%             
%             [ ypos(i), j ] = max( [ ypos1, ypos2 ] );
%             if j == 1
%                 ipos(i) = ipos1;
%             else
%                 ipos(i) = ipos2;
%             end;
%             ypos(i) = ypos(i)*10;
% %             ypos = ylim(2)/10;
%             
%         else
%             
%             if length( mx_plot{i,1} ) > 2;
%                 ipos1 = ceil( ( length(mx_plot{i,1}) + 1 ) / 2 );
%                 ypos1 = mx_plot{i,1}(ipos1);
%             else
%                 [ ypos1, ipos1 ] = min( mx_plot{i,1} );
%             end
%             
%             if length( mx_plot{i,2} ) > 2
%                 ipos2 = ceil( ( length(mx_plot{i,2}) + 1 ) / 2 );
%                 ypos2 = mx_plot{i,2}(ipos2);
%             else
%                 [ ypos2, ipos2 ] = min( mx_plot{i,2} );
%             end
%             
%             [ ypos(i), j ] = min( [ ypos1, ypos2 ] );
%             if j == 1
%                 ipos(i) = ipos1;
%             else
%                 ipos(i) = ipos2;
%             end;
%             ypos(i) = ypos(i)/10;
% %             ypos(i) = ylim(2)/100;
%             
%         end
%        
%         if j == 1;
%             xpos(i) = aa_plot{i,1}(ipos(i));
%         else
%             xpos(i) = aa_plot{i,2}(ipos(i));
%         end
%
%         margin = 2 * ( pos(4)/saxis_pos(iplot,4) );
%         ypos(i) = min( ylim(2)/margin, max( ylim(1)*margin, ypos(i) ) );
% 
%         ipos(i) = find( aa_plot{i} == 2*iz );
%         if isempty(ipos(i))
%             ipos(i) = 1;
%         end
%         
%         xpos(i) = aa_plot{i}(ipos(i));

        if strcmp(yscale,'log')
            if tmp_flag
                ypos(i) = max( ypos( xpos == xpos(i) ) ) * 3;
            else
                ypos(i) = min( ypos( xpos == xpos(i) ) ) / 3;
            end
            ypos(i) = min( ylim(2)/3, max( ylim(1)*3, ypos(i) ) );
        else
            dypos = diff(ylim)/5;
            if tmp_flag
                ypos(i) = max( ypos( xpos == xpos(i) ) ) + dypos;
            else
                ypos(i) = min( ypos( xpos == xpos(i) ) ) - dypos;
            end
            ypos(i) = min( ylim(2) - dypos, max( ylim(1) + dypos, ypos(i) ) );
        end
        
        for iplot = 1:nsubplot
            
            if plot_flag{i}(iplot)
                
                if nmx == 3
                    mask = mx_plot{i,1} >= mx_min;
                end
                    
                jj = 1;
                for j = 1:nmx
                   
                    if nmx ~= 3
                        mask = mx_plot{i,j} >= mx_min;
                    end
                    
                    if any( ~mask )
                        
                        plot( saxis_h(iplot), aa_plot{i}(~mask), mx_plot{i,j}(~mask), ...
                                'DisplayName', element{iz+1}, ...
                                'LineStyle', 'none', ...
                                'LineWidth', LineWidth(jj), ...
                                'Marker', MarkerSet{jj,:}, ...
                                'MarkerSize', MarkerSize(jj), ...
                                'MarkerFaceColor', 'none', ...
                                'Color', ColorSet(k,:), ...
                                'Parent', group_h(iplot,1) );
                        
                        if ~isempty( aa_plot{i}(mask) )
                            mx_h{j,iplot}(i) = plot( saxis_h(iplot), aa_plot{i}(mask), mx_plot{i,j}(mask), ...
                                'DisplayName', element{iz+1}, ...
                                'LineStyle', LineSet{jj,:}, ...
                                'LineWidth', LineWidth(jj), ...
                                'Marker', MarkerSet{jj,:}, ...
                                'MarkerSize', MarkerSize(jj), ...
                                'MarkerFaceColor', 'none', ...
                                'Color', ColorSet(k,:), ...
                                'Parent', group_h(iplot,2) );
                        end
                    else
                        mx_h{j,iplot}(i) = plot( saxis_h(iplot), aa_plot{i}, mx_plot{i,j}, ...
                            'DisplayName', element{iz+1}, ...
                            'LineStyle', LineSet{jj,:}, ...
                            'LineWidth', LineWidth(jj), ...
                            'Marker', MarkerSet{jj,:}, ...
                            'MarkerSize', MarkerSize(jj), ...
                            'MarkerFaceColor', 'none', ...
                            'Color', ColorSet(k,:), ...
                            'Parent', group_h(iplot,2) );
                    end
                    
                    jj = mod( jj, nlines ) + 1;
                    
                end
                
                text_h{iplot}(i) = text( double(xpos(i)), double(ypos(i)), element{iz+1}, ...
                    'Parent', saxis_h(iplot), ...
                    'Color', ColorSet(k,:), ...
                    'FontSize', 18, ...
                    'Interpreter', 'latex', ...
                    'Clipping', 'on', ...
                    'HorizontalAlignment', 'center',...
                    'VerticalAlignment', valign );
                
            end 
        end
        
        k = mod( k, ncolors ) + 1;
        tmp_flag = ~tmp_flag;
        
    end
    
end

varargout{1} = zz_plot;
varargout{2} = aa_plot;
varargout{3} = mx_plot;
varargout{4} = saxis_h;
varargout{5} = text_h;
varargout{6} = mx_h;