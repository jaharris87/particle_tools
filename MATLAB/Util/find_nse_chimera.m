function [ time_nse, nse_flag, nse_index ] = find_nse_chimera( time, radius, theta, time_nse_chimera, radius_nse, theta_nse, plist )
% Determine when particles are in NSE by using NSE radius from frame files
    
    % Initialize
    time_nse    = zeros( size(plist) );
    nse_flag    = true( size(radius) );
    nse_index   = ones( size(plist) );
    ny          = length( theta_nse );
    
    % Define edges of theta grid
    theta_edges = linspace( eps, pi, ny+1 );
    
    % Only need inner edges
    theta_edges = theta_edges(1:end-1);
    
    if length(time_nse) > 1
        
        % Find the frame corresponding to the particle data
        [~,i_time] = histc( time, [time_nse_chimera;Inf] );
        i_time(i_time == 0) = 1;
        
        % Determine to which ray the particles belongs
        [~,i_theta_bin] = histc( theta, theta_edges );
        
        for i = 1:length(time)
            
            % Find the frame corresponding to the particle data at time(i)
%             i_nse = find( time_nse_chimera <= time(i), 1, 'last' );
%             if isempty(i_nse)
%                 i_nse = 1;
%             end
%             if i_time(i) ~= i_nse
%                 warning('%s%d','histc does not match find for i_nse: ',i);
%             end
            i_nse = i_time(i);
            
            % Determine to which ray the particle belongs
%             [~,i_theta] = histc( theta(i,plist), theta_edges );
%             if any( i_theta ~= i_theta_bin(i,plist) )
%                 warning('%s%d','histc does not match find for i_theta: ',i);
%             end
            i_theta = i_theta_bin(i,plist);
            
            for j = 1:ny
                
                p_bin = plist(i_theta==j);
                
                % Compare the radius to the NSE radius to see if the particle is in NSE
                nse_flag(i,p_bin(radius(i,p_bin) > radius_nse(i_nse,j))) = false;
                
                % Find the last instance in time of the particle being in NSE
                nse_index(p_bin(radius(i,p_bin) <= radius_nse(i_nse,j))) = i;
                
            end
            
            % Handle glitched particles separately
            p_glitch             = plist(i_theta==ny+1 | i_theta==0);
            nse_flag(i,p_glitch) = true;
            nse_index(p_glitch)  = length(time);
            
        end
        
        time_nse(plist(nse_index > 0)) = time(nse_index(plist(nse_index > 0)));
        
    else
        
        % Find the frame corresponding to the particle data
        [~,i_time] = histc( time, [time_nse_chimera;Inf] );
        i_time(i_time == 0) = 1;
        
        % Determine to which ray the particles belongs
        [~,i_theta_bin] = histc( theta, theta_edges );      
        
        for i = 1:length(time)
            
            % Find the frame corresponding to the particle data at time(i)
%             i_nse = find( time_nse_chimera <= time(i), 1, 'last' );
%             if isempty(i_nse)
%                 i_nse = 1;
%             end
%             if i_time(i) ~= i_nse
%                 warning('%s%d','histc does not match find for i_nse: ',i);
%             end
            i_nse = i_time(i);
            
            % Determine to which ray the particle belongs
%             [~,i_theta] = histc( theta(i), theta_edges );
%             if i_theta ~= i_theta_bin(i)
%                 warning('%s%d','histc does not match find for i_theta: ',i);
%             end
            i_theta = i_theta_bin(i);
            
            % Compare the radius to the NSE radius to see if the particle is in NSE
            if i_theta == 0
                nse_flag(i) = true;
                nse_index   = length(time);
            elseif radius(i) > radius_nse(i_nse,i_theta)
                nse_flag(i) = false;
            else
                nse_index = i;
            end
            
        end
        
        if nse_index > 0
            time_nse = time(nse_index);
        end
        
    end

end