function [ x_plot, y_plot, out_mask ] = downsample_plot( change_min )

fig_h  = copyobj(gcf,0);

% copyobj(get(fig_ho,'children'),fig_h);
axis_h = findobj(fig_h,'Type','axes');
for k = 1:length(axis_h)
    if strcmp(get(axis_h(k),'YScale'),'log')
        logplot = true;
    else
        logplot = false;
    end
    line_h = findobj(axis_h(k),'Type','line');
    if ~isempty(line_h)
%         obj_h = line_h;
        for m = 1:length(line_h);
            obj_h = line_h(m);
            xdata = get(obj_h(1),'XData');
            ydata = cell(size(obj_h));
            for j = 1:length(obj_h)
                ydata{j} = get(obj_h(j),'YData');
            end
            ntime = length(xdata);
            
            if logplot
                for j = 1:length(obj_h)
                    ydata{j} = log10(ydata{j});
                end
            end
            
            change    = zeros(size(obj_h));
            r_old     = zeros(size(obj_h));
            maxchange = zeros(size(xdata));
            
            for i = 1:ntime
                for j = 1:length(obj_h)
                    change(j) = abs( 1.0 - ydata{j}(i)*r_old(j) );
                end
                maxchange(i) = max(change);
                if maxchange(i) >= change_min || i == 1 || i == ntime
                    for j = 1:length(obj_h)
                        if isnan(ydata{j}(i))
                            r_old(j) = Inf;
                        else
                            r_old(j) = 1/ydata{j}(i);
                        end
                    end
                end
            end
            out_mask        = ( maxchange >= change_min );
            out_mask(1)     = true;
            out_mask(ntime) = true;
            
            out_mask(1:end-1) = out_mask(1:end-1) | out_mask(2:end);
            
            x_plot = xdata(out_mask);
            disp(length(x_plot));
            for j = 1:length(obj_h)
                if logplot
                    y_plot = 10.^ydata{j}(out_mask);
                else
                    y_plot = ydata{j}(out_mask);
                end
                
                set(obj_h(j),'XData',x_plot,'YData',y_plot);
            end
        end
    end
end

end

