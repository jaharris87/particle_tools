function [ time_nse, nse_flag, nse_index ] = find_nse_index(time,temp,plist,temp_nse)
    time_nse    = zeros( size(plist) );
    nse_flag    = true( size(temp) );
    nse_index   = ones( size(plist) );
    temp_peak   = max(temp);
    for p_id = plist
        if temp_peak(p_id) > temp_nse && temp(end,p_id) <= temp_nse
            nse_index(p_id) = find( temp(:,p_id) > temp_nse, 1, 'last' );
            nse_flag(nse_index(p_id)+1:length(time),p_id) = false;
        elseif temp(end,p_id) >= temp_nse
            nse_index(p_id) = size(temp,1);            
        elseif temp_peak(p_id) < temp_nse
            nse_index(p_id)  = 1;
            nse_flag(:,p_id) = false;
        end
    end
    time_nse(plist(nse_index > 0)) = time(nse_index(plist(nse_index > 0)));
end