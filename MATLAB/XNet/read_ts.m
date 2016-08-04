function [ pp_xn, pp_xntot, pp_time, pp_temp, pp_rho, pp_zz, pp_aa, i_alpha, pp_xninitial, pp_xnfinal, pp_flx ] = read_ts( ts_filename_base, plist, tracer_mass, time )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

i = 1;
ts_filename = sprintf('%s%04.0f', ts_filename_base, i);

% Call Raph's ts-file reader
[pp_zz, pp_aa, ~, ~, ~, ~, ~, ~ ] = read_ts_file( ts_filename );

i_he4 = find(pp_aa==4 & pp_zz==2);
i_ni56 = find(pp_aa==56 & pp_zz==28);
i_si28 = find(pp_aa==28 & pp_zz==14);
i_fe56 = find(pp_aa==56 & pp_zz==26);
i_fe52 = find(pp_aa==52 & pp_zz==26);
i_cr48 = find(pp_aa==48 & pp_zz==24);
i_ti44 = find(pp_aa==44 & pp_zz==22);
i_s32 = find(pp_aa==32 & pp_zz==16);
i_c12 = find(pp_aa==12 & pp_zz==6);
i_o16 = find(pp_aa==16 & pp_zz==8);
i_ne20 = find(pp_aa==20 & pp_zz==10);
i_mg24 = find(pp_aa==24 & pp_zz==12);
i_ar36 = find(pp_aa==36 & pp_zz==18);
i_ca40 = find(pp_aa==40 & pp_zz==20);
i_zn60 = find(pp_aa==60 & pp_zz==30);
i_nn = find(pp_aa==1 & pp_zz==0);
i_pp = find(pp_aa==1 & pp_zz==1);

i_alpha = [i_he4,i_c12,i_o16,i_ne20,i_mg24,i_si28,i_s32,i_ar36,i_ca40,i_ti44,i_cr48,i_fe52,i_ni56,i_zn60];

% Interpolate ts-file mass fractions to match particle times
pp_xn = cell(1,length(plist));
pp_xninitial = zeros(length(pp_zz),length(plist));
pp_xnfinal = zeros(length(pp_zz),length(plist));
pp_xntot = zeros(max(pp_aa),length(plist));
pp_flx = cell(1,length(plist));
pp_time = cell(1,length(plist));
pp_temp = cell(1,length(plist));
pp_rho = cell(1,length(plist));

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
    
%     pp_xn = zeros(length(time),length(pp_zz));
    pp_xntot_tmp = zeros(1,max(pp_aa));
    
    ts_filename = sprintf('%s%04.0f', ts_filename_base, p_id);
%     ts_filename = sprintf('%s%02.0f', ts_filename_base, id);
    disp(ts_filename);
    
    [~, ~, pp_xn{id}, pp_time{id}, pp_temp{id}, pp_rho{id}, ~, ~ ] = read_ts_file( ts_filename );
    
    [~,ia,~]    = unique(pp_time{id});
    pp_time{id} = pp_time{id}(ia);
    pp_temp{id} = pp_temp{id}(ia);
    pp_rho{id}  = pp_rho{id}(ia);
    pp_xn{id}   = pp_xn{id}(:,ia);
%     pp_flx{id}  = tmp_flx;
    
    pp_xninitial(:,id) = pp_xn{id}(:,1);
    pp_xnfinal(:,id)   = pp_xn{id}(:,end);
    
    for i = 1:max(pp_aa)
        pp_xntot_tmp(i) = sum( pp_xn{id}(pp_aa==i,end) );
    end
    
    pp_xntot(:,id) = pp_xntot_tmp;
    
%     if length(pp_time{id}) > 1
%         for j = 1:length(pp_zz)
%             pp_x(:,j) = interp1( pp_time{id}, pp_xn(j,:), time, 'linear', pp_xn(j,end) );
%         end
%     else
%         pp_x = repmat(pp_xn',length(time),1);
%     end
%     pp_xntot = pp_xntot + pp_x*tracer_mass(p_id);
%     pp_xninitial(:,p_id) = pp_xn(1,:);
%     pp_xnfinal(:,p_id) = pp_xn(end,:);
        
%     % Update progress monitor
%     ppm.increment(id);
    
end
% for p_id = plist

% % Clean-up for progress monitor
% delete(ppm);
%     
% delete(poolobj);

% pp_xn = zeros(length(time),length(pp_zz));

% for i = 1:length(p_nonnse_reduced_theta)
%     ts_filename = sprintf('%s%03.0f','C:\Users\Austin\Documents\research\particles\B12-WH07\results\alpha-np-fe56_reduced_theta\ts\ts_',i);
%     disp(ts_filename);
%     [~, ~, pp_x, pp_time, ~, ~, ~, ~, ~, ~ ] = read_ts_file( ts_filename );
% %     xn = zeros(length(time),length(pp_zz));
%     xn = interp1( pp_time, pp_x', time, 'linear', 0.0 );
%     xn_sum_reduced_theta = xn_sum_reduced_theta + xn*tracer_mass(i);
% end


% for i = 1:length(p_nonnse_extra)
%     for j = 1:length(pp_zz)
%         xn(:,j,p_nonnse_extra(i)) = interp1( pp_time{i}, pp_x{i}(j,:), time, 'linear', 0.0 );
%     end
% end
% xn_sum_SN150_8GK_noextrap=sum(xn(:,:,p_nonnse_extra(1:end-40))*M_tracer,3);
% xn_sum_SN150_8GK_noextrap=xn_sum_SN150_8GK_noextrap + sum(xn(:,:,p_nonnse_extra(end-39:end))*M_tracer*0.400807528404510,3);

% save('B12-WH07\pp_data.mat', 'pp_*','xn','xn_sum');


end

