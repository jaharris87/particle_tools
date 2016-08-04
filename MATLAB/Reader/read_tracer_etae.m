function [ etae,notfound ] = read_tracer_etae( etae_fname_base, plist, time_input )

time_start = time_input(1);
time_end = [0,0];

num_fields = 2;
notfound(plist) = 0.0;

% Determine time range of data files and adjust time range to fit data if necessary
for p_id = 1:length(plist)
    
    time_min = time_start;
    time_max = time_input(end);
    
    etae_fname = sprintf('%s%05.0f', etae_fname_base, p_id);
    etae_fform = repmat('%f',1,num_fields);
    
    etae_id = fopen(etae_fname);
    if etae_id >= 0
        file_read = textscan(etae_id, etae_fform, 1);
        file_read = cell2mat(file_read);
        time_init = file_read(1, 1);
        
        offset = -ftell(etae_id) - 1;
        fseek(etae_id, offset, 'eof');
        file_read = textscan(etae_id, etae_fform, 1);
        file_read = cell2mat(file_read);
        time_final = file_read(1, 1);
        
        ii = 1;
        jj = length(time_input);
        while time_min < time_init && ii < jj
            ii = ii + 1;
            time_min = time_input(ii);
        end
        
        while time_max > time_final && jj > ii
            jj = jj - 1;
            time_max = time_input(jj);
        end
        
        if time_min > time_start
            time_start = time_min;
        end
        if time_final > time_end(1)
            time_end = [time_max, time_final];
        end
        
        fclose(etae_id);
    else
        notfound(p_id) = p_id;
    end
      
end

ii = find(time_input==time_start);
jj = find(time_input==time_end(1));

% time_input = [time_input(ii:jj-1),time_end];
% time_input = [time_input(ii:jj-1),time_max];
time_input = time_input(ii:end);

etae     = zeros(length(time_input),length(plist));

poolobj = parpool;
tic;

parfor p_id = plist
    
    time_range = sort(time_input);
    
    % Generate filename from etae_fname_base and partilce ID
    etae_fname = sprintf('%s%05.0f', etae_fname_base, p_id);
    disp(etae_fname);

    % Open the file
    etae_id = fopen(etae_fname);
    if etae_id >= 0
        
        % Read from large file in chunks until end-of-file or stop time
        etae_fform = repmat('%f',1,num_fields);
        etae_data  = textscan(etae_id, etae_fform);
        etae_data  = cell2mat(etae_data);
        
        % Eliminate duplicate lines
        [~,ia,~] = unique( etae_data(:,1) );
        etae_data  = etae_data(ia,:,:);
        
        if( length(time_range) > 1)
            % Interpolate data points to common points in time.
            p_data = interp1( etae_data(:,1), etae_data, time_range, 'linear', 'extrap' );
        else
            p_data = repmat(etae_data(1,:,:),[length(time_range),1,1]);
        end
        
        % Write temporary particle data in array for plotting later
        etae     (:,p_id) = p_data(:,2);
        
        % Close the file
        fclose(etae_id);
    end
    
    etae_data=[];
end
elapsed_time=toc;
display(elapsed_time);

delete(poolobj);

end