function [ time, nu_e0, nu_e1, nu_e2, nu_e3, nu_eta ] = read_tracer_nu( nu_fname_base, plist, varargin )

num_fields = 21;

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
        
        nu_fname = sprintf('%s%05.0f', nu_fname_base, p_id);
        nu_fform = repmat('%f',1,num_fields);
        
        nu_id = fopen(nu_fname);
        if nu_id >= 0
            file_read = textscan(nu_id, nu_fform, 1);
            file_read = cell2mat(file_read);
            time_init = file_read(1, 1);
            
            offset = -ftell(nu_id) - 1;
            fseek(nu_id, offset, 'eof');
            file_read = textscan(nu_id, nu_fform, 1);
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
            
            fclose(nu_id);
        end
        
    end
    
    ii = find(time_interp==time_start);
    jj = find(time_interp==time_end(1));
    
    % time_interp = [time_interp(ii:jj-1),time_end];
    % time_interp = [time_interp(ii:jj-1),time_max];
    time_interp = unique( [ time_init, time_interp(ii:end) ] );
    
    time   = time_interp';
    nu_e0  = zeros(length(time_interp),4,length(plist));
    nu_e1  = zeros(length(time_interp),4,length(plist));
    nu_e2  = zeros(length(time_interp),4,length(plist));
    nu_e3  = zeros(length(time_interp),4,length(plist));
    nu_eta = zeros(length(time_interp),4,length(plist));
    
    poolobj = parpool;
    tic;
    
    parfor p_id = 1:length(plist)
        
        time_range = sort(time_interp);
        
        % Generate filename from nu_fname_base and partilce ID
        nu_fname = sprintf('%s%05.0f', nu_fname_base, p_id);
        disp(nu_fname);
        
        % Open the file
        nu_id = fopen(nu_fname);
        if nu_id >= 0
            
            % Read from large file in chunks until end-of-file or stop time
            nu_fform = repmat('%f',1,num_fields);
            nu_data  = textscan(nu_id, nu_fform);
            nu_data  = cell2mat(nu_data);
            
            % Eliminate duplicate lines
            [~,ia,~] = unique( nu_data(:,1) );
            nu_data  = nu_data(ia,:);
            
            if( length(time_range) > 1)
                % Interpolate data points to common points in time.
                p_data = interp1( nu_data(:,1), nu_data, time_range, 'linear', 'extrap' );
            else
                p_data = repmat(nu_data(1,:,:),[length(time_range),1,1]);
            end
            
            % Write temporary particle data in array for plotting later
            for i = 1:4
                nu_e0 (:,i,p_id) = p_data(:,1+i);
                nu_e1 (:,i,p_id) = p_data(:,5+i);
                nu_e2 (:,i,p_id) = p_data(:,9+i);
                nu_e3 (:,i,p_id) = p_data(:,13+i);
                nu_eta(:,i,p_id) = p_data(:,17+i);
            end
            
            % Close the file
            fclose(nu_id);
        end
        
        nu_data=[];
    end
    elapsed_time=toc;
    display(elapsed_time);
    
    delete(poolobj);
    
else
    
    p_id = plist(1);
    
    % Generate filename from nu_fname_base and partilce ID
    nu_fname = sprintf('%s%05.0f', nu_fname_base, p_id);
    
    % Open the file
    nu_id = fopen(nu_fname);
    if nu_id >= 0
        
        % Read from large file in chunks until end-of-file or stop time
        nu_fform = repmat('%f',1,num_fields);
        nu_data  = textscan(nu_id, nu_fform);
        nu_data  = cell2mat(nu_data);
        
        % Eliminate duplicate lines
        [~,ia,~] = unique( nu_data(:,1) );
        nu_data  = nu_data(ia,:);
        
        % Write temporary particle data in array for plotting later
        time   = nu_data(:,1);
        nu_e0  = zeros(length(time),4);
        nu_e1  = zeros(length(time),4);
        nu_e2  = zeros(length(time),4);
        nu_e3  = zeros(length(time),4);
        nu_eta = zeros(length(time),4);
        for i = 1:4
            nu_e0 (:,i) = nu_data(:,1+i);
        end
        for i = 1:4
            nu_e1 (:,i) = nu_data(:,5+i);
        end
        for i = 1:4
            nu_e2 (:,i) = nu_data(:,9+i);
        end
        for i = 1:4
            nu_e3 (:,i) = nu_data(:,13+i);
        end
        for i = 1:4
            nu_eta(:,i) = nu_data(:,17+i);
        end
        
        % Close the file
        fclose(nu_id);
    end
end

nu_eta(nu_eta<0.0) = 0.0;

end