function [ time, radius, theta, v_rad, v_theta, temp, density, ye, enpy, pe_int, pe_bind, press, lapse, dpe_nuc, dpe_neut ] = read_tracer_fate( fname_base, plist, varargin )

num_fields = 15;

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
    for p_id = plist
        
        time_min = time_start;
        time_max = time_interp(end);
        
        fate_fname = sprintf('%s%05.0f', fname_base, p_id);
        fate_fform = repmat('%f',1,num_fields);
        
        f_id = fopen(fate_fname);
        file_read = textscan(f_id, fate_fform, 1);
        file_read = cell2mat(file_read);
        time_init = file_read(1, 1);
        
        offset = -ftell(f_id) - 1;
        fseek(f_id, offset, 'eof');
        file_read = textscan(f_id, fate_fform, 1);
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
        
        fclose(f_id);
        
    end

    ii = find(time_interp==time_start);
    jj = find(time_interp==time_end(1));
%     time_interp = [time_interp(ii:jj-1),time_end];
%     time_interp = [time_interp(ii:jj-1),time_max];
    time_interp = unique( [ time_init, time_interp(ii:end) ] );
%     display(time_interp(1:2));
%     display(time_interp(end-1:end));
    
    time     = time_interp';
    radius   = zeros(length(time_interp),length(plist));
    theta    = zeros(length(time_interp),length(plist));
    v_rad    = zeros(length(time_interp),length(plist));
    v_theta  = zeros(length(time_interp),length(plist));
    temp     = zeros(length(time_interp),length(plist));
    density  = zeros(length(time_interp),length(plist));
    ye       = zeros(length(time_interp),length(plist));
    enpy     = zeros(length(time_interp),length(plist));
    pe_int   = zeros(length(time_interp),length(plist));
    pe_bind  = zeros(length(time_interp),length(plist));
    press    = zeros(length(time_interp),length(plist));
    lapse    = zeros(length(time_interp),length(plist));
    dpe_nuc  = zeros(length(time_interp),length(plist));
    dpe_neut = zeros(length(time_interp),length(plist));
    
    poolobj = parpool;
    tic;
    
    parfor p_id = 1:length(plist)
        
        time_range=sort(time_interp);
        fate_fform = repmat('%f',1,num_fields);
        
        % Generate filename from fname_base and partilce ID
        fate_fname = sprintf('%s%05.0f', fname_base, p_id);
        disp(fate_fname);
        
        % Open the file
        f_id = fopen(fate_fname);
        
        % Read from large file in chunks until end-of-file or stop time
        f_data = textscan( f_id, fate_fform );
        f_data = cell2mat( f_data );
        
        % Eliminate duplicate lines
        [~,ia,~] = unique( f_data(:,1) );
        f_data = f_data(ia,:);
        
        if( length(f_data(:,1)) > 1)
            % Interpolate data points to common points in time.
            p_data = interp1( f_data(:,1), f_data, time_range, 'linear', 'extrap' );
        else
            p_data = repmat(f_data(1,:),length(time_range),1);
        end
        
        % Write temporary particle data in array for plotting later
        radius  (:,p_id) = p_data(:,2);
        theta   (:,p_id) = p_data(:,3);
        v_rad   (:,p_id) = p_data(:,4);
        v_theta (:, p_id) = p_data(:,5);
        temp    (:, p_id) = p_data(:,6);
        density (:, p_id) = p_data(:,7);
        ye      (:, p_id) = p_data(:,8);
        enpy    (:, p_id) = p_data(:,9);
        pe_int  (:, p_id) = p_data(:,10);
        pe_bind (:, p_id) = p_data(:,11);
        press   (:, p_id) = p_data(:,12);
        lapse   (:, p_id) = p_data(:,13);
        dpe_nuc (:, p_id) = p_data(:,14);
        dpe_neut (:,p_id) = p_data(:,15);
        
        % Close the file
        fclose(f_id);
        
        f_data=[];
    end
    
    elapsed_time=toc;
    display(elapsed_time);
    delete(poolobj);

else
    
    p_id = plist(1);
    
    fate_fname = sprintf('%s%05.0f', fname_base, p_id);
    fate_fform = repmat('%f',1,num_fields);
    
    f_id = fopen(fate_fname);
    
    % Read from large file in chunks until end-of-file or stop time
    f_data = textscan( f_id, fate_fform );
    f_data = cell2mat( f_data );
    
    % Eliminate duplicate lines
    [~,ia,~] = unique( f_data(:,1) );
    f_data = f_data(ia,:);
    
    time      = f_data(:,1);
    radius    = f_data(:,2);
    theta     = f_data(:,3);
    v_rad     = f_data(:,4);
    v_theta   = f_data(:,5);
    temp      = f_data(:,6);
    density   = f_data(:,7);
    ye        = f_data(:,8);
    enpy      = f_data(:,9);
    pe_int    = f_data(:,10);
    pe_bind   = f_data(:,11);
    press     = f_data(:,12);
    lapse     = f_data(:,13);
    dpe_nuc   = f_data(:,14);
    dpe_neut  = f_data(:,15);
    
    fclose(f_id);
    
end

c = 2.99792458E+10;  % Speed of light [cm s^{-1}]

radius(radius<0.0)   = eps;
theta(theta>pi)      = pi;
theta(theta<0.0)     = 0.0;
temp(temp<0.0)       = eps;
v_rad(v_rad>=c)      = eps;
v_theta(v_theta>=c)  = eps;
density(density<0.0) = eps;
press(press<0.0)     = eps;
ye(ye<0.0)           = eps;
ye(ye>1)             = 1.0;

end