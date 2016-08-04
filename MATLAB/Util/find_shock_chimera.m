function [ time_shock, shock_flag, shock_index ] = find_shock_chimera( time, radius, theta, time_shock_chimera, radius_shock, theta_shock, plist )
% Determine when particles are in NSE by using NSE radius from frame files
    
    % Initialize
    time_shock    = zeros( 1, size(radius,2) );
    shock_flag    = true( size(radius) );
    shock_index   = uint16( zeros( 1, size(radius,2) ) );
    ny            = length( theta_shock );
    
    % Define edges of theta grid
    theta_edges = linspace( eps, pi, ny+1 );
    
    % Only need inner edges
    theta_edges = theta_edges(1:end-1);
    
    if length(time_shock) > 1
        
        % Find the frame corresponding to the particle data
        [~,i_time] = histc( time, [time_shock_chimera;Inf] );
        i_time(i_time == 0) = 1;
        
        % Determine to which ray the particles belongs
        [~,i_theta_bin] = histc( theta, theta_edges );
        
        for i = length(time):-1:1
            
            % Find the frame corresponding to the particle data at time(i)
%             i_shock = find( time_shock_chimera <= time(i), 1, 'last' );
%             if isempty(i_shock)
%                 i_shock = 1;
%             end
%             if i_time(i) ~= i_shock
%                 warning('%s%d','histc does not match find for i_shock: ',i);
%             end
            i_shock = i_time(i);
            
            % Determine to which ray the particle belongs
%             [~,i_theta] = histc( theta(i,plist), theta_edges );
%             if any( i_theta ~= i_theta_bin(i,plist) )
%                 warning('%s%d','histc does not match find for i_theta: ',i);
%             end
            i_theta = i_theta_bin(i,plist);
            
            for j = 1:ny
                
                p_bin = plist(i_theta==j);
                
                % Compare the radius to the NSE radius to see if the particle is in NSE
                shock_flag(i,p_bin(radius(i,p_bin) > radius_shock(i_shock,j))) = false;
                
                % Find the last instance in time of the particle being in NSE
                shock_index(p_bin(radius(i,p_bin) <= radius_shock(i_shock,j))) = i;
                
            end
            
            % Handle glitched particles separately
            p_glitch             = plist(i_theta==ny+1 | i_theta==0);
            shock_flag(i,p_glitch) = true;
            shock_index(p_glitch)  = length(time);
            
        end
        
        for id = 1:length(plist)
            p_id = plist(id);
            if shock_index(p_id) > 0
                shock_flag(shock_index(p_id):end,p_id) = true;
            end
        end
        
        time_shock(plist(shock_index > 0)) = time(shock_index(plist(shock_index > 0)));
        
    else
        
        % Find the frame corresponding to the particle data
        [~,i_time] = histc( time, [time_shock_chimera;Inf] );
        i_time(i_time == 0) = 1;
        
        % Determine to which ray the particles belongs
        [~,i_theta_bin] = histc( theta, theta_edges );      
        
        for i = length(time):-1:1
            
            % Find the frame corresponding to the particle data at time(i)
%             i_shock = find( time_shock_chimera <= time(i), 1, 'last' );
%             if isempty(i_shock)
%                 i_shock = 1;
%             end
%             if i_time(i) ~= i_shock
%                 warning('%s%d','histc does not match find for i_shock: ',i);
%             end
            i_shock = i_time(i);
            
            % Determine to which ray the particle belongs
%             [~,i_theta] = histc( theta(i), theta_edges );
%             if i_theta ~= i_theta_bin(i)
%                 warning('%s%d','histc does not match find for i_theta: ',i);
%             end
            i_theta = i_theta_bin(i);
            
            % Compare the radius to the NSE radius to see if the particle is in NSE
            if i_theta == 0
                shock_flag(i) = true;
                shock_index   = length(time);
            elseif radius(i) > radius_shock(i_shock,i_theta)
                shock_flag(i) = false;
            else
                shock_index = i;
            end
            
        end
        
        if shock_index > 0
            shock_flag(shock_index:end) = true;
            time_shock = time(shock_index);
        end
        
    end

end