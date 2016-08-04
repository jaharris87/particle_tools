function [ time_nse, nse_flag, nse_index ] = find_nse_temp( time, density, temp, plist )
% Determine when particles are in NSE by using T_NSE equation from B-series
    
    % Initialize
    time_nse    = zeros( 1, size(density,2) );
    nse_flag    = true( size(density) );
    nse_index   = uint16( zeros( 1, size(density,2) ) );
    
    
    rho_h = 2.0e8;
    rho_l = 5.0e7;
    t_h   = 6.5e9;
    t_l   = 5.7e9;
    slope = ( t_h - t_l ) / ( rho_h - rho_l );
    
    temp_nse = ( density - rho_l ) * slope + t_l;
    temp_nse(temp_nse>t_h) = t_h;
    
    % Convert to GK
    temp_nse  = temp_nse * 1.0e-9;
    temp_dnse = temp_nse - 0.2;
    
    if length(time_nse) > 1

        for i = 1:length(plist)
            
            p_id = plist(i);
            
            % See if the particle ever reaches NSE
            if all( temp(:,p_id) < temp_nse(:,p_id) )
                nse_flag(1:length(time),p_id) = false;
                nse_index(p_id) = 0;
                
                % Otherwise, see if the particle ever drops out of NSE
            elseif all( temp(:,p_id) > temp_dnse(:,p_id) )
                nse_flag(1:length(time),p_id) = true;
                nse_index(p_id) = length(time);
            else
                
                % Find the last time the particle "deflashes"
                i_nse = find( temp(:,p_id) > temp_nse(:,p_id), 1, 'last' );
                i_deflash = find( temp(i_nse:end,p_id) < temp_dnse(i_nse:end,p_id), 1, 'first' );
                
                % Make sure the particle deflashed after the final flash
                if ~isempty(i_deflash)
                    i_deflash = i_deflash + i_nse - 1;
                    nse_flag(i_deflash:length(time),p_id) = false;
                    nse_index(p_id) = i_deflash - 1;
                    
                    % Otherwise, consider the particle as always in NSE from a post-processing perspective
                else
                    nse_flag(1:length(time),p_id) = true;
                    nse_index(p_id) = length(time);
                end
            end
        end
        
        time_nse(plist(nse_index > 0)) = time(nse_index(plist(nse_index > 0)));
        
    else
        
        % See if the particle ever reaches NSE
        if all( temp(:) < temp_nse(:) )
            nse_flag(1:length(time)) = false;
            nse_index = 0;
            
            % Otherwise, see if the particle ever drops out of NSE
        elseif all( temp(:) > temp_dnse(:) )
            nse_flag(1:length(time)) = true;
            nse_index = length(time);
        else
            
            % Find the last time the particle "deflashes"
            i_nse = find( temp(:) > temp_nse(:), 1, 'last' );
            i_deflash = find( temp(i_nse:end) < temp_dnse(i_nse:end), 1, 'first' );
            
            % Make sure the particle deflashed after the final flash
            if ~isempty(i_deflash)
                i_deflash = i_deflash + i_nse - 1;
                nse_flag(i_deflash:length(time)) = false;
                nse_index = i_deflash - 1;
                
                % Otherwise, consider the particle as always in NSE from a post-processing perspective
            else
                nse_flag(1:length(time)) = true;
                nse_index = length(time);
            end
        end
        
        if nse_index > 0
            time_nse = time(nse_index);
        end
        
    end
    
end