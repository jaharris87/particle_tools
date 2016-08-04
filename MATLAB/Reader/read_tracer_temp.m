function [ time, temp, density, ye, flxtot, nu_temp ] = read_tracer_temp( temp_fname_base, plist, varargin )

num_fields = 12;

if nargin > 2
    time_interp = varargin{1};
    interp_flag = true;
else
    interp_flag = false;
end

if interp_flag

    time_start = time_interp(1);
    time_end = [0,0];
    
    % Determine time range of data files and adjust time range to fit data if necessary
    for p_id = 1:length(plist)
        
        time_min = time_start;
        time_max = time_interp(end);
        
        temp_fname = sprintf('%s%05.0f', temp_fname_base, p_id);
        temp_fform = repmat('%f',1,num_fields);
        
        temp_id = fopen(temp_fname);
        if temp_id >= 0
            file_read = textscan(temp_id, temp_fform, 1);
            file_read = cell2mat(file_read);
            time_init = file_read(1, 1);
            
            offset = -ftell(temp_id) - 1;
            fseek(temp_id, offset, 'eof');
            file_read = textscan(temp_id, temp_fform, 1);
            file_read = cell2mat(file_read);
            time_final = file_read(1, 1);
            
            ii = 1;
            jj = length(time_interp);
            while time_min < time_init && ii < jj
                ii = ii + 1;
                time_min = time_interp(ii);
            end
            
            while time_max > time_final && jj > ii
                jj = jj - 1;
                time_max = time_interp(jj);
            end
            
            if time_min > time_start
                time_start = time_min;
            end
            if time_final > time_end(1)
                time_end = [time_max, time_final];
            end
            
            fclose(temp_id);
        end
        
    end
    
    ii = find(time_interp==time_start);
    jj = find(time_interp==time_end(1));
    
    % time_interp = [time_interp(ii:jj-1),time_end];
    % time_interp = [time_interp(ii:jj-1),time_max];
    time_interp = unique( [ time_init, time_interp(ii:end) ] );
    
    time     = time_interp';
    temp     = zeros(length(time_interp),length(plist));
    density  = zeros(length(time_interp),length(plist));
    ye       = zeros(length(time_interp),length(plist));
    flxtot   = zeros(length(time_interp),4,length(plist));
    nu_temp  = zeros(length(time_interp),4,length(plist));
    
    poolobj = parpool;
    tic;
    
    parfor p_id = 1:length(plist)
        
        time_range = sort(time_interp);
        
        % Generate filename from temp_fname_base and partilce ID
        temp_fname = sprintf('%s%05.0f', temp_fname_base, p_id);
        disp(temp_fname);
        
        % Open the file
        temp_id = fopen(temp_fname);
        if temp_id >= 0
            
            % Read from large file in chunks until end-of-file or stop time
            temp_fform = repmat('%f',1,num_fields);
            temp_data  = textscan(temp_id, temp_fform);
            temp_data  = cell2mat(temp_data);
            
            % Eliminate duplicate lines
            [~,ia,~] = unique( temp_data(:,1) );
            temp_data  = temp_data(ia,:);
            
            if( length(time_range) > 1)
                % Interpolate data points to common points in time.
                p_data = interp1( temp_data(:,1), temp_data, time_range, 'linear', 'extrap' );
            else
                p_data = repmat(temp_data(1,:,:),[length(time_range),1,1]);
            end
            
            % Write temporary particle data in array for plotting later
            temp     (:,p_id) = p_data(:,2);
            density  (:,p_id) = p_data(:,3);
            ye       (:,p_id) = p_data(:,4);
            for i = 1:4
                flxtot   (:,i,p_id) = p_data(:,4+i);
                nu_temp  (:,i,p_id) = p_data(:,8+i);
            end
            
            % Close the file
            fclose(temp_id);
        end
        
        temp_data=[];
    end
    elapsed_time=toc;
    display(elapsed_time);
    
    delete(poolobj);
    
else
    
    p_id = plist(1);
    
    % Generate filename from temp_fname_base and partilce ID
    temp_fname = sprintf('%s%05.0f', temp_fname_base, p_id);
    
    % Open the file
    temp_id = fopen(temp_fname);
    if temp_id >= 0
        
        % Read from large file in chunks until end-of-file or stop time
        temp_fform = repmat('%f',1,num_fields);
        temp_data  = textscan(temp_id, temp_fform);
        temp_data  = cell2mat(temp_data);
        
        % Eliminate duplicate lines
        [~,ia,~] = unique( temp_data(:,1) );
        temp_data  = temp_data(ia,:);
        
        % Write temporary particle data in array for plotting later
        time     = temp_data(:,1);
        temp     = temp_data(:,2);
        density  = temp_data(:,3);
        ye       = temp_data(:,4);
        
        flxtot       = zeros(length(time),4);
        nu_temp      = zeros(length(time),4);
        for i = 1:4
            flxtot   (:,i) = temp_data(:,4+i);
            nu_temp  (:,i) = temp_data(:,8+i);
        end
        
        % Close the file
        fclose(temp_id);
    end
end

temp(temp<0.0)       = eps;
density(density<0.0) = eps;
ye(ye<0.0)           = eps;
ye(ye>1)             = 1.0;
nu_temp(nu_temp<0.0) = eps;
nu_temp(isnan(nu_temp)) = eps;

end