function [ nu_flx, nu_lum, nu_rms ] = read_tracer_flux( flux_fname_base, fate_fname_base, plist, time_input, radius_input )

unui( 1) = 2.5749E+00;
unui( 2) = 5.0118E+00;
unui( 3) = 6.2304E+00;
unui( 4) = 7.7452E+00;
unui( 5) = 9.6284E+00;
unui( 6) = 1.1969E+01;
unui( 7) = 1.4880E+01;
unui( 8) = 1.8498E+01;
unui( 9) = 2.2995E+01;
unui(10) = 2.8586E+01;
unui(11) = 3.5537E+01;
unui(12) = 4.4177E+01;
unui(13) = 5.4918E+01;
unui(14) = 6.8271E+01;
unui(15) = 8.4870E+01;
unui(16) = 1.0551E+02;
unui(17) = 1.3116E+02;
unui(18) = 1.6305E+02;
unui(19) = 2.0269E+02;
unui(20) = 2.5197E+02;

ergmev  = 1.6021773E-06;   % Ergs per MeV

num_fields = 4;
% A = cell(numel(plist));

nu_flx   = zeros(length(time_input),4,length(plist));
nu_lum   = zeros(length(time_input),4,length(plist));
nu_rms   = zeros(length(time_input),4,length(plist));

% matlabpool
tic;

for p_id = plist
    
    time_range = sort(time_input);

    fate_fname = sprintf('%s%05.0f', fate_fname_base, p_id);
    fate_id    = fopen(fate_fname);
    fate_fform = repmat('%f',1,25);
    fate_data  = textscan(fate_id, fate_fform);
    fate_data  = cell2mat(fate_data);
    [~,ia,~]   = unique( fate_data(:,1) );
    fate_data  = fate_data(ia,:);
    
    % Generate filename from fname_base and partilce ID
    flux_fname = sprintf('%s%05.0f', flux_fname_base, p_id);
    display(flux_fname);

    % Open the file
    flux_id = fopen(flux_fname);

    % Read from large file in chunks until end-of-file or stop time
    flux_fform = repmat('%f',1,4);
    flux_data = textscan( flux_id, flux_fform );
    flux_data = cell2mat( flux_data );
    flux_data = reshape( flux_data', [], 4, 20 );
    
    % Eliminate duplicate lines
%     [~,ia,~] = unique( flux_data(:,1) );
    flux_data = flux_data(ia,:,:);
    
    if( length(time_range) > 1)
        % Interpolate data points to common points in time.
        p_data = interp1( fate_data(:,1), flux_data, time_range, 'linear', 'extrap' );
    else
        p_data = repmat(flux_data(1,:,:),[length(time_range),1,1]);
    end
    
    % Write temporary particle data in array for plotting later
    for i = 1:4
        for k = 1:20
            nu_flx(:,i,p_id) = nu_flx(:,i,p_id) + p_data(:,i,k);
            nu_lum(:,i,p_id) = nu_lum(:,i,p_id) + p_data(:,i,k)*unui(k); 
            nu_rms(:,i,p_id) = nu_rms(:,i,p_id) + p_data(:,i,k)*unui(k)^2;
        end
        nu_rms(:,i,p_id) = sqrt(nu_rms(:,i,p_id) ./ nu_flx(:,i,p_id));
        nu_lum(:,i,p_id) = 4 * pi * radius_input(:,p_id).^2 * ergmev * 1.0e-51;
    end
    
    % Close the file
    fclose(flux_id);
    fclose(fate_id);
    
    flux_data=[];
end
elapsed_time=toc;
display(elapsed_time);
% matlabpool close

end