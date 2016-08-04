function [ time_nse, radius_nse, theta_nse ] = read_nse_radius( nse_fname_base, ny, start_frame, end_frame )

num_fields = 3;
    
time_nse   = zeros(end_frame,1);
radius_nse = zeros(end_frame,ny);

theta_nse = linspace( pi / (2*ny), (2*ny - 1)*pi / (2*ny), ny );

poolobj = parpool;
tic;

parfor frame = start_frame:end_frame
    
    % Generate filename from temp_fname_base and partilce ID
    nse_fname = sprintf('%s%05.0f.dat', nse_fname_base, frame);
    disp(nse_fname);
    
    % Open the file
    nse_id = fopen(nse_fname);
    if nse_id >= 0
        
        % Read the time for this frame, skipping header lines
        time_fform = repmat('%f',1,2);
        time_data  = textscan(nse_id, time_fform, 1, 'HeaderLines', 5);
        time_data  = cell2mat(time_data);
        
        time_nse(frame) = time_data(:,1);
        
        % Read the NSE radius for each theta zone
        for j = 1:3
            fgets(nse_id);
        end
        
        nse_fform = strcat('%d',repmat('%f',1,num_fields));
        nse_data  = textscan(nse_id, nse_fform);
        nse_data  = nse_data(2:end);
        nse_data  = cell2mat(nse_data);
        
        % Write temporary particle data in array for plotting later
        radius_nse(frame,:) = nse_data(:,3);
        
        % Close the file
        fclose(nse_id);
    end
    
    time_data=[];
    nse_data=[];
end
elapsed_time=toc;
display(elapsed_time);

delete(poolobj);

end