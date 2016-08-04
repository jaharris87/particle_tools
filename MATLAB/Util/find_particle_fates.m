function [ p_pns, p_notpns, p_shocked, p_unshocked, p_bound, p_unbound, p_posvr, p_negvr, p_nse_chimera, p_nonnse_chimera, p_nse_8GK, p_nonnse_8GK, varargout ] = find_particle_fates( time, radius, temp, density, pe_kin, pe_thermal, pe_grav, v_rad, nse_flag_chimera, shock_index, plist )
% Determine fates of particles

    % Determine "boundedness" ratio
    e_grav = -( pe_kin + pe_thermal ) ./ pe_grav;
    
    % Determine when each particle hits the shock
%     [ ~, shock_flag, shock_index ] = find_shock_chimera( time, radius, theta, chim_time_shock, chim_radius_shock, chim_theta_shock, plist );
    
    % Determine when each particle is in NSE or non-NSE
%     [ ~, nse_flag_chimera, ~ ] = find_nse_chimera( time, radius, theta, chim_time_nse, chim_radius_nse, chim_theta_nse, plist );
    [ ~, nse_flag_8GK, ~ ] = find_nse_index( time, temp, plist, 8 );
%     [ ~, nse_flag_chimera, ~ ] = find_nse_temp( time, density, temp, plist );

    % Create the cells which will contain particle IDs of each categorization
    p_pns               = cell(length(time),1);
    p_notpns            = cell(length(time),1);
    p_shocked           = cell(length(time),1);
    p_unshocked         = cell(length(time),1);
    p_bound             = cell(length(time),1);
    p_unbound           = cell(length(time),1);
    p_posvr             = cell(length(time),1);
    p_negvr             = cell(length(time),1);
    p_nse_chimera       = cell(length(time),1);
    p_nonnse_chimera    = cell(length(time),1);
    p_nse_8GK           = cell(length(time),1);
    p_nonnse_8GK        = cell(length(time),1);
    
    p_data              = zeros(size(pe_grav));
    
    for i = 1:length(time)
        
        % Determine which particles are in/out the PNS
        p_pns{i}    = uint16(plist( or(density(i,:) >= 1e11, radius(i,:) < 3e6) ...
                                 | and(density(i,:) >= 1e10, radius(i,:) < 6e6) ));
        p_notpns{i} =  uint16(plist( setdiff(plist,p_pns{i}) ));
        
        % Of the particles not in PNS, determine which of them have been shocked/not shocked
        p_shocked{i} = uint16(p_notpns{i}(shock_index(p_notpns{i}) > 0 & shock_index(p_notpns{i}) <= i));
        if ~isempty(p_shocked{i}) 
            % Of the partices which have reached the shock, determine which of them are bound/unbound
            p_bound{i}   = uint16(p_shocked{i}(e_grav(i,p_shocked{i}) <= 1 & shock_index(p_shocked{i}) <= i-200 ));
            p_unbound{i} = uint16(p_shocked{i}(e_grav(i,p_shocked{i}) > 1 ));
            
            if ~isempty(p_unbound{i})
                % Of the particles which are unbound, determine which are outwardly/inwardly streaming
                p_posvr{i} = uint16(p_unbound{i}( v_rad(i,p_unbound{i}) > 0 ));
                p_negvr{i} = uint16(p_unbound{i}( v_rad(i,p_unbound{i}) < 0 & shock_index(p_unbound{i}) <= i-200));
            end
        end
        
        % In addition to particles which have not yet been shocked, include particles that haven't been shocked long enough to judge boundedness
        p_unshocked{i} = uint16(plist( setdiff(p_notpns{i}, unique( [p_bound{i},p_unbound{i},p_posvr{i},p_negvr{i}]) )));
        
        p_nse_chimera{i}    = uint16(p_notpns{i}(nse_flag_chimera(i,p_notpns{i})));
        p_nonnse_chimera{i} = uint16(p_notpns{i}(~nse_flag_chimera(i,p_notpns{i})));
        p_nse_8GK{i}        = uint16(p_notpns{i}(nse_flag_8GK(i,p_notpns{i})));
        p_nonnse_8GK{i}     = uint16(p_notpns{i}(~nse_flag_8GK(i,p_notpns{i})));
        
        p_data(i,p_bound{i})     = 1;
        p_data(i,p_negvr{i})     = 2;
        p_data(i,p_posvr{i})     = 3;
        p_data(i,p_unshocked{i}) = 4;
        p_data(i,p_pns{i})       = 5;
        
    end
    
    varargout{1} = p_data;
    
end