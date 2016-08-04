function [ time_shock, radius_shock, theta_shock ] = read_shock_radius( shock_fname_base, ny, start_frame, end_frame )

num_fields = 2;
    
time_shock   = zeros(end_frame,1);
radius_shock = zeros(end_frame,ny);

theta_shock = linspace( pi / (2*ny), (2*ny - 1)*pi / (2*ny), ny );

poolobj = parpool;
tic;

parfor frame = start_frame:end_frame
    
    % Generate filename
    shock_fname = sprintf('%s%05.0f.dat', shock_fname_base, frame);
    disp(shock_fname);
    
    % Open the file
    shock_id = fopen(shock_fname);
    if shock_id >= 0
        
        % Read the time for this frame, skipping header lines
        time_fform = '# %f %f';
        time_data  = textscan(shock_id, time_fform, 1, 'HeaderLines', 5);
        time_data  = cell2mat(time_data);
        
        time_shock(frame) = time_data(:,1);
        
        % Read the shock radius for each theta zone
        for j = 1:2
            fgets(shock_id);
        end
        
        shock_fform = repmat('%f',1,num_fields);
        shock_data  = textscan(shock_id, shock_fform);
        shock_data  = cell2mat(shock_data);
        
        radius_shock(frame,:) = shock_data(:,1);
        
        % Close the file
        fclose(shock_id);
    end
    
    time_data=[];
    shock_data=[];
end
elapsed_time=toc;
display(elapsed_time);

delete(poolobj);

end