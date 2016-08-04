function [ th_time, th_temp, th_rho, th_ye, th_nueflux, th_nuebarflux, th_numutauflux, th_numutaubarflux ] = read_th_file( th_fname )
% function [ th_time, th_temp, th_rho, th_ye, th_nueflux, th_nuebarflux ] = read_th_file( th_fname_base, plist )
% function [ th_time, th_temp, th_rho, th_ye ] = read_th_file( th_fname_base, plist )

th_fform = '%d %f %f %f %f %f %f %f %f %f %f %f %f';

% for p_id = plist
%     th_fname = sprintf( '%s%05.0f', th_fname_base, p_id );
%     disp(th_fname);
    f_id = fopen( th_fname );
    textscan( f_id, '%s' , 5, 'Delimiter', '\n' );
    th_data = textscan( f_id, th_fform );
    th_time = th_data{2};
    th_temp = th_data{3};
    th_rho = th_data{4};
    th_ye = th_data{13};
    th_nueflux = th_data{5}*1e-42;
    th_nuebarflux = th_data{6}*1e-42;
    th_numutauflux = th_data{7}*1e-42;
    th_numutaubarflux = th_data{8}*1e-42;
    fclose( f_id );
    clear th_data;
% end