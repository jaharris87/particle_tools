function [ shock_index, shocked, varargout ] = find_shocked_particles( time, time_bounce, press, density, pe_int, pe_kin, plist )
% Determine when particles become shocked by finding discontinuity in mach number

    % Array index corresponding to bounce
    [~, bounce_index] = min( abs( time_bounce-time ) );

    % Adiabatic exponent
    gamma = 1 ./ (1 - press./(density.*pe_int));
    
    % Sound speed
    v_sound(:,plist) = sqrt( gamma(:,plist) .* press(:,plist) ./ density(:,plist) );
    
    % Mach number
    mach = sqrt( 2*pe_kin ) ./ v_sound;

    % Initialize
    shock_index(plist) = uint16(length(time)+1);
    shocked = false(size(press));

    % Find discontinuity in Mach number
    %   Treat the Mach number as a "signal" and use analysis by wavelet decomposition to detect shock.
    %   Haar wavelets are used due to the fact that short wavelets are often more effective than long ones
    %   at detecting signal ruptures.
    for p_id = plist
        
        result = [];
        
        % Compute the multilevel 1-D wavelet decomposition at level 1
        [C,L] = wavedec(mach(:,p_id),1,'haar');
        
        % Extract the detail coefficients at level 1 from the wavelet decomposition
        D = detcoef(C,L,1);
        
        % Interpolate across the time domain
        D = abs(interpft(D,length(time)));
        maxpeak = max(D(bounce_index:end));
        
        % Find peaks in the detail coefficients, corresponding to discontinuities
        if( maxpeak > 0.005 )
            
            % Isolate the large peaks and choose the first peak as the shock discontinuity
            [~,locs] = findpeaks(D(bounce_index:end),'sortstr','descend','minpeakheight',0.005,'minpeakdistance',200);
            
            result = min(locs);
            
        end
        
        % Store the array index of the shock discontinuity for each particle
        if ~isempty(result)
            shock_index(p_id) = uint16(result+bounce_index-1);
            shocked(shock_index(p_id):numel(mach(:,1)),p_id) = 1;
        end
        
    end
    
    % Optional output from intermediate results
    nout = max(nargout,2) - 2;
    switch( nout )
        case 1
            varargout{1} = bounce_index;
        case 2
            varargout{1} = bounce_index;
            varargout{2} = gamma;
        case 3
            varargout{1} = bounce_index;
            varargout{2} = gamma;
            varargout{3} = v_sound;
        case 4
            varargout{1} = bounce_index;
            varargout{2} = gamma;
            varargout{3} = v_sound;
            varargout{4} = mach;
    end

end