function [ pp_zz, pp_aa, pp_xn, pp_time, pp_temp, pp_rho, pp_xn0, pp_xnf, varargout ] = read_ts_files( ts_filename_base, plist )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ndigits = length(sprintf('%d',max(plist)));
% ndigits = 1;
filename_format = sprintf( '%s%d%s','%s%0',ndigits,'.0f' );

i = 1;
ts_filename = sprintf(filename_format, ts_filename_base, i);

% Call Raph's ts-file reader
[pp_zz, pp_aa, ~, ~, ~, ~, ~, ~ ] = read_ts_file( ts_filename );

% Interpolate ts-file mass fractions to match particle times
pp_time = cell(1,length(plist));
pp_temp = cell(1,length(plist));
pp_rho = cell(1,length(plist));
pp_xn = cell(1,length(plist));
pp_xn0 = zeros(length(pp_zz),length(plist));
pp_xnf = zeros(length(pp_zz),length(plist));
pp_flx = cell(1,length(plist));

% poolobj = parpool;

% % Initialize progress monitor
% N = length(plist);
% percentage_update = 0.05;
% do_debug = 0;
% run_javaaddpath = 1;
% show_execution_time = 1;
% ppm = ParforProgressStarter2('read_ts', N, percentage_update, do_debug, run_javaaddpath, show_execution_time);

% parfor id = 1:length(plist)
for id = 1:length(plist)
    
    p_id = plist(id);
    
    ts_filename = sprintf(filename_format, ts_filename_base, p_id);
    
    disp(ts_filename);
    
    [~, ~, pp_xn{id}, pp_time{id}, pp_temp{id}, pp_rho{id}, ~, ~, ~ ] = read_ts_file( ts_filename );
    
    [~,ia,~]    = unique(pp_time{id});
    pp_time{id} = pp_time{id}(ia);
    pp_temp{id} = pp_temp{id}(ia);
    pp_rho{id}  = pp_rho{id}(ia);
    pp_xn{id}   = pp_xn{id}(:,ia);
    if nargout > 8
        pp_flx{id} = tmp_flx;
    end
    
    pp_xn0(:,id) = pp_xn{id}(:,1);
    pp_xnf(:,id) = pp_xn{id}(:,end);
        
%     % Update progress monitor
%     ppm.increment(id);
    
end

% % Clean-up for progress monitor
% delete(ppm);
%     
% delete(poolobj);

if nargout > 8
    varargout{1} = pp_flx;
end

% for i = 1:length(p_nonnse_extra)
%     for j = 1:length(pp_zz)
%         xn(:,j,p_nonnse_extra(i)) = interp1( pp_time{i}, pp_x{i}(j,:), time, 'linear', 0.0 );
%     end
% end
% xn_sum_SN150_8GK_noextrap=sum(xn(:,:,p_nonnse_extra(1:end-40))*M_tracer,3);
% xn_sum_SN150_8GK_noextrap=xn_sum_SN150_8GK_noextrap + sum(xn(:,:,p_nonnse_extra(end-39:end))*M_tracer*0.400807528404510,3);

% save('B12-WH07\pp_data.mat', 'pp_*','xn','xn_sum');


end

