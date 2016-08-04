%% Initialization/Setup
if ispc
    homedir = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
else
    homedir = getenv('HOME');
end

model_mass = '12';
model_name = strcat('B',model_mass,'-WH07');
load_constants_B12;

% Create lists of particle IDs and any extras necessary to cover the grid in Chimera
plist = 1:4000;
p_grid = uint16(reshape(plist,40,[]));
p_extra = 4001:4720;
plist_all = [ plist, p_extra ];
tracer_mass(plist_all) = M_tracer;
% Adjust the mass of the last row if particles overlap the grid boundary
tracer_mass(plist_all(end-39:end)) = M_tracer * (4.13753599987209e+33 - 4.13158E+33) / (4.14644E+33 - 4.13158E+33);

% WH07 progenitor properties
wh07_presn_nspecies   = 1612;
wh07_presn_start_line = 2;
wh07_presn_end_line   = 1100;
wh07_nuc_start_line   = 2;
wh07_nuc_end_line     = 764;
wh07_fin_start_line   = 2;
wh07_fin_end_line     = 764;

wh07_presn_offgrid   = 451:1098;
wh07_offgrid         = 115:762;

% Relative path for tracer_reader ASCII data
fate_fname_base  = fullfile( model_name, 'full_profile', strcat('fate_S',model_mass,'-') );
nse_fname_base   = fullfile( model_name, 'nse_radius',   'nse_radius' );
shock_fname_base = fullfile( model_name, 'shock_p',      'shock_p' );
temp_fname_base  = fullfile( model_name, 'full_profile', strcat('temp_S',model_mass,'-') );
nu_fname_base    = fullfile( model_name, 'full_profile', strcat('nu_S',model_mass,'-') );

% Relative path for progenitor data
prad_rows_fname      = fullfile( model_name, 'pinput_file' );
wh07_species_fname   = fullfile( model_name, 'raw_species_final.d' );
wh07data_presn_fname = fullfile( homedir, 'chimera', 'trunk', 'Initial_Models', '2NSE-2nonNSE_Pre_Collapse_Model_Generator', 'progenitors', 'wh07',     strcat('s',model_mass,'_wh07_presn.d') );
wh07data_nuc_fname   = fullfile( homedir, 'chimera', 'trunk', 'Initial_Models', '2NSE-2nonNSE_Pre_Collapse_Model_Generator', 'progenitors', 'wh07_fin', strcat('s',model_mass,'_wh07_nucleo.d') );
wh07data_fin_fname   = fullfile( homedir, 'chimera', 'trunk', 'Initial_Models', '2NSE-2nonNSE_Pre_Collapse_Model_Generator', 'progenitors', 'wh07_fin', strcat('s',model_mass,'_wh07_final.d') );

wh07data_nucwind_fname = fullfile( homedir, 'chimera', 'trunk', 'Initial_Models', '2NSE-2nonNSE_Pre_Collapse_Model_Generator', 'progenitors', 'wh07_fin', strcat('s',model_mass,'_wh07_nucleo_wind.d') );
wh07data_finwind_fname = fullfile( homedir, 'chimera', 'trunk', 'Initial_Models', '2NSE-2nonNSE_Pre_Collapse_Model_Generator', 'progenitors', 'wh07_fin', strcat('s',model_mass,'_wh07_final_wind.d') );

% Create MAT-files for MATLab data
model_fname      = strcat( model_name, '.mat' );
model_work_fname = strcat( model_name, '_work.mat' );

raw_fname        = fullfile( model_name, 'raw_data.mat' );
energy_fname     = fullfile( model_name, 'energy_data.mat' );
nse_fname        = fullfile( model_name, 'nse_data.mat' );
shock_fname      = fullfile( model_name, 'shock_data.mat' );
pdata_fname      = fullfile( model_name, 'pfate_data.mat' );
mass_fname       = fullfile( model_name, 'mass_data.mat' );
nu_fname         = fullfile( model_name, 'neutrino_data.mat' );
pp_fname         = fullfile( model_name, 'postprocessing_data.mat' );
chim_fname       = fullfile( model_name, 'chimera_data.mat' );
grid_fname       = fullfile( model_name, 'grid_data.mat' );
wh07_fname       = fullfile( model_name, 'wh07_data.mat' );
wh07_presn_fname = fullfile( model_name, 'wh07_presn_data.mat' );
wh07_nuc_fname   = fullfile( model_name, 'wh07_nuc_data.mat' );
wh07_fin_fname   = fullfile( model_name, 'wh07_fin_data.mat' );
solar_fname      = 'solar.mat';

raw        = matfile( raw_fname,        'Writable', true );
pe         = matfile( energy_fname,     'Writable', true );
nse        = matfile( nse_fname,        'Writable', true );
shock      = matfile( shock_fname,      'Writable', true );
p          = matfile( pdata_fname,      'Writable', true );
m          = matfile( mass_fname,       'Writable', true );
nu         = matfile( nu_fname,         'Writable', true );
pp         = matfile( pp_fname,         'Writable', true );
chim       = matfile( chim_fname,       'Writable', true );
grid       = matfile( grid_fname,       'Writable', true );
wh07       = matfile( wh07_fname,       'Writable', true );
wh07_presn = matfile( wh07_presn_fname, 'Writable', true );
wh07_nuc   = matfile( wh07_nuc_fname,   'Writable', true );
wh07_fin   = matfile( wh07_fin_fname,   'Writable', true );
solar      = matfile( solar_fname,      'Writable', true );

%% Raw particle data from Chimera simulation

% Create a reduced sample of time domain to be shared amongst all particles for analysis
time_input = get_time_input( 0.0, time_final, time_bounce, 5000, 10 );

[ time, radius, theta, v_rad, v_theta, temp, density, ye, enpy, pe_int, pe_bind, press, lapse, dpe_nuc, dpe_neut ] ...
    = read_tracer_fate( fate_fname_base, plist, time_input );

% Get some important time indices
[~,itime_bounce]    = min( abs( time_bounce-time ) );
[~,itime]           = min(abs(time-time_bseries));
[~,itime_m50]       = min(abs(time-(time(itime)-0.05)));
[~,itime_m75]       = min(abs(time-(time(itime)-0.075)));
[~,itime_m100]      = min(abs(time-(time(itime)-0.1)));
[~,itime_m150]      = min(abs(time-(time(itime)-0.15)));
[~,itime_500]       = min(abs(time-0.5-time_bounce));
[~,itime_600]       = min(abs(time-0.6-time_bounce));
[~,itime_700]       = min(abs(time-0.7-time_bounce));
[~,itime_800]       = min(abs(time-0.8-time_bounce));
[~,itime_900]       = min(abs(time-0.9-time_bounce));
[~,itime_1000]      = min(abs(time-1.0-time_bounce));

% Identify particles that have left the grid and adjust physical quantities to last known on-grid value
itime_offgrid(plist_all) = length(time);
for p_id = plist
    if any(radius(:,p_id) > max_radius)
        itime_offgrid(p_id) = find( radius(:,p_id) <= max_radius, 1, 'last' ) + 1;
        
        theta(itime_offgrid(p_id):end,p_id)     = theta(itime_offgrid(p_id)-1,p_id);
        v_rad(itime_offgrid(p_id):end,p_id)     = v_rad(itime_offgrid(p_id)-1,p_id);
        v_theta(itime_offgrid(p_id):end,p_id)   = v_theta(itime_offgrid(p_id)-1,p_id);
        temp(itime_offgrid(p_id):end,p_id)      = temp(itime_offgrid(p_id)-1,p_id);
        density(itime_offgrid(p_id):end,p_id)   = density(itime_offgrid(p_id)-1,p_id);
        pe_int(itime_offgrid(p_id):end,p_id)    = pe_int(itime_offgrid(p_id)-1,p_id);
        pe_bind(itime_offgrid(p_id):end,p_id)   = pe_bind(itime_offgrid(p_id)-1,p_id);
        press(itime_offgrid(p_id):end,p_id)     = press(itime_offgrid(p_id)-1,p_id);
        lapse(itime_offgrid(p_id):end,p_id)     = lapse(itime_offgrid(p_id)-1,p_id);
        dpe_nuc(itime_offgrid(p_id):end,p_id)   = dpe_nuc(itime_offgrid(p_id)-1,p_id);
        dpe_neut(itime_offgrid(p_id):end,p_id)  = dpe_neut(itime_offgrid(p_id)-1,p_id);
        radius(itime_offgrid(p_id):end,p_id)    = radius(itime_offgrid(p_id)-1,p_id);
    end
end

save( raw_fname, ...
    'time','radius','theta','v_rad','v_theta','density','temp','ye','press','enpy','lapse','dpe_neut','dpe_nuc' );

%% Energy data
M_r        = find_enclosed_mass( radius, M_inner, M_tracer );
pe_grav    = - G * M_r ./ radius;
pe_kin     = 0.5 * (v_rad.^2 + v_theta.^2);
pe_thermal = pe_int - pe_int_offset - me*ye*ku - pe_bind + ye*dmnp;

save( energy_fname, ...
    'pe_int','pe_grav','pe_kin','pe_thermal','pe_bind' );

%% NSE data
[ chim_time_nse, chim_radius_nse, chim_theta_nse ] = read_nse_radius( nse_fname_base, ny_frame, start_frame, end_frame );
[ time_nse_chimera, nse_flag_chimera, itime_nse_chimera ] = find_nse_chimera( time, radius, theta, chim_time_nse, chim_radius_nse, chim_theta_nse, plist );

[ time_nse_8GK, nse_flag_8GK, itime_nse_8GK ] = find_nse_index( time, temp, plist, 8 );

save( nse_fname, ...
    'chim_time_nse','chim_radius_nse','chim_theta_nse','time_nse_*','nse_flag_*','itime_nse_*' );

%% Shock data
[ chim_time_shock, chim_radius_shock, chim_theta_shock ] = read_shock_radius( shock_fname_base, ny_frame, start_frame, end_frame );

[~,iframe_bseries] = min(abs(chim_time_shock-time_bseries));
iframe_offgrid(1:ny_frame) = end_frame;
chim_dt = [0;diff(chim_time_shock(:))];
for itheta = 1:ny_frame
    chim_dr = [chim_radius_shock(1,itheta); diff(chim_radius_shock(:,itheta))];
    chim_vr = chim_dr ./ chim_dt;
    rtest = [ chim_radius_shock(1,itheta); chim_radius_shock(1:end-1,itheta) + chim_vr(1:end-1) .* chim_dt(2:end) ];
    if any( rtest > max_radius );
        iframe_offgrid(itheta) = find( rtest > max_radius, 1, 'first' ) - 1;
        chim_radius_shock(iframe_offgrid(itheta):end,itheta) = max_radius;
    end
end
clear chim_dt chim_dr chim_vr rtest

[ time_shock, shock_flag, shock_index ] = find_shock_chimera( time, radius, theta, chim_time_shock, chim_radius_shock, chim_theta_shock, plist );

save( shock_fname, ...
    'chim_time_shock','chim_radius_shock','chim_theta_shock','time_shock','shock_flag','shock_index' );

%% Particle fate data
[p_pns, p_notpns, p_shocked, p_unshocked, p_bound, p_unbound, p_posvr, p_negvr, p_nse_chimera, p_nonnse_chimera, p_nse_8GK, p_nonnse_8GK, p_data] ...
    = find_particle_fates( time, radius, temp, density, ...
                           pe_kin, pe_thermal, pe_grav, v_rad, ...
                           nse_flag_chimera, shock_index, plist );
for i = 1:length(time)
    p_unshocked{i} = [ p_unshocked{i}, p_extra ];
end
p_data(:,p_extra) = 4;

save( pdata_fname, ...
    'p_*', 'tracer_mass' );

%% Progenitor grid data
[ me_grid, rade_grid, vr_grid, rho_grid, temp_grid, press_grid, e_grid, enpy_grid, vth_grid, abar_grid, ye_grid ] = import_grid_wh07( wh07data_presn_fname );
me_grid   = [ 0.0; me_grid ];
rade_grid = [ 0.0; rade_grid ];
radc_grid = 0.5 * ( rade_grid(1:end-1) + rade_grid(2:end) );

i_gridedge = find( rade_grid <= max_radius, 1, 'last' );

save( grid_fname, ...
    '*_grid','i_gridedge' );

%% Enclosed mass data
M_pmasscut = sum( tracer_mass( [p_pns{end},p_bound{end}] ) ) + M_inner;
i_masscut  = find( me_grid <= M_masscut, 1, 'last' );
i_pmasscut = find( me_grid <= M_pmasscut, 1, 'last' );

r_cut  = rade_grid(i_masscut);
r_pcut = rade_grid(i_pmasscut);

for i = 1:length(time)
    [~,i_m100] = min(abs(time-(time(i)-0.1)));
    m_pns(i) = sum(tracer_mass(p_pns{i}))/M_solar; 
    m_unb(i) = sum(tracer_mass(p_unbound{i}))/M_solar;
    m_bound(i) = sum(tracer_mass(p_bound{i}))/M_solar;
    m_posvr(i) = sum(tracer_mass(p_posvr{i}))/M_solar;
    m_negvr(i) = sum(tracer_mass(p_negvr{i}))/M_solar;
    m_negvr100(i) = sum(tracer_mass(intersect( unique( [p_negvr{i_m100:i}] ),p_unbound{i})))/M_solar;
    m_unshocked(i) = sum(tracer_mass(p_unshocked{i}))/M_solar;
end

save( mass_fname, ...
    'm_*','M_*','i_masscut','i_pmasscut','r_cut','r_pcut' );

%% Neutrino data
[ ~, ~, ~, ~, flxtot, nu_temp ] = read_tracer_temp( temp_fname_base, plist, time_input );
for p_id = plist(itime_offgrid<length(time))
    for i = 1:4
        flxtot(itime_offgrid(p_id):end,i,p_id) = flxtot(itime_offgrid(p_id)-1,i,p_id) .* ( radius(itime_offgrid(p_id)-1,p_id) ./ radius(itime_offgrid(p_id):end,p_id) ).^2;
        nu_temp(itime_offgrid(p_id):end,i,p_id) = nu_temp(itime_offgrid(p_id)-1,i,p_id);
    end
end

save( nu_fname, ...
    'flxtot','nu_temp' );

%% Chimera simulation data
chim_time = chim_time_nse;
[~,chim_itime] = min(abs(chim_time-time_bseries));

save( chim_fname, ...
    'chim_time','chim_itime' );

%% Solar abundance data
nname_solar = solar.nname_solar;
aa_solar = solar.aa_solar;
zz_solar = solar.zz_solar;
xn_solar = solar.xn_solar;
y_solar  = solar.y_solar;
inuc_za_solar = solar.inuc_za_solar;
inuc_a_solar = solar.inuc_a_solar;
inuc_z_solar = solar.inuc_z_solar;

xn_a_solar = zeros( 1, max(aa_solar) );
xn_z_solar = zeros( 1, max(zz_solar) );
y_a_solar = zeros( 1, max(aa_solar) );
y_z_solar = zeros( 1, max(zz_solar) );

for i = 1:max(aa_solar)
    xn_a_solar(i) = sum( xn_solar(inuc_a_solar(i)) );
    y_a_solar (i) = sum( y_solar (inuc_a_solar(i)) );
end

for i = 1:max(zz_solar)
    xn_z_solar(i) = sum( xn_solar(inuc_z_solar(i)) );
    y_z_solar (i) = sum( y_solar (inuc_z_solar(i)) );
end

%% Post-processing data

% Functions for relative composition differences 
delta  = @(x1,x2) log10(x1./x2);
rdelta = @(x1,x2) sum(abs(log10(x1(x1>0 & x2>0)./x2(x1>0 & x2>0)))) ./ sum(abs(0.5*log10(x1(x1>0 & x2>0).*x2(x1>0 & x2>0))));

save( pp_fname, ...
    'delta','rdelta','tracer_mass' );

%% WH07 data
[ wh07_nname, wh07_zz, wh07_aa ] = import_wh07_species( wh07_species_fname );
wh07_nspecies = length(wh07_zz);
[ solar_in_wh07, wh07_solar, wh07_in_solar, solar_wh07 ] = map_nuclei( wh07_aa, wh07_zz, aa_solar, zz_solar );

% Pre-SN data
[ wh07_presn_m, wh07_presn_dm, wh07_presn_r, ~, wh07_presn_vr, wh07_presn_rho, ~, wh07_presn_temp, wh07_presn_press, wh07_presn_eint, wh07_presn_enpy, wh07_presn_vth, wh07_presn_abar, wh07_presn_ye ,wh07_presn_xn, ~, ~, ~, ~ ] ...
    = import_wh07( wh07data_presn_fname, wh07_presn_nspecies, wh07_presn_start_line, wh07_presn_end_line );

wh07_presn_mx = wh07_presn_xn .* repmat( wh07_presn_dm, 1, wh07_presn_nspecies );

save( wh07_presn_fname, ...
    'wh07_presn_*' );

% 100 second data
wh07_nuc_nspecies = wh07_nspecies;
[ wh07_m, wh07_dm, wh07_nuc_r, ~, wh07_nuc_vr, wh07_nuc_rho, ~, wh07_nuc_temp, wh07_nuc_press, wh07_nuc_eint, wh07_nuc_enpy, wh07_nuc_vth, wh07_nuc_abar, wh07_nuc_ye, wh07_nuc_xn, ~, ~, ~, ~ ] ...
    = import_wh07( wh07data_nuc_fname, wh07_nuc_nspecies, wh07_nuc_start_line, wh07_nuc_end_line );

inuc_z_wh07  = @(zz) get_iso( wh07_aa, wh07_zz, 'Z', zz );
inuc_a_wh07  = @(aa) get_iso( wh07_aa, wh07_zz, 'A', aa );
inuc_za_wh07 = @(zz,aa) get_iso( wh07_aa, wh07_zz, 'A', aa, 'Z', zz );

save( wh07_fname, ...
    'wh07_aa','wh07_zz','wh07_solar','wh07_nname','wh07_nspecies','wh07_offgrid','wh07_m','wh07_dm','solar_wh07','wh07_in_solar','solar_in_wh07','inuc_*_wh07' );

wh07_nuc_offgrid  = wh07_offgrid;
wh07_nuc_m  = wh07_m;
wh07_nuc_dm = wh07_dm;
wh07_nuc_mx = wh07_nuc_xn .* repmat( wh07_nuc_dm, 1, wh07_nuc_nspecies );

wh07_nuc_mx_wind = import_wh07_wind( wh07data_nucwind_fname, wh07_nuc_nspecies );
wh07_nuc_mxtotal_ejecta  = sum( wh07_nuc_mx(2:end,:), 1 ) + wh07_nuc_mx_wind;
wh07_nuc_xntotal_ejecta  = wh07_nuc_mxtotal_ejecta ./ sum( wh07_nuc_mxtotal_ejecta );
wh07_nuc_ytotal_ejecta   = wh07_nuc_xntotal_ejecta ./ wh07_aa';
wh07_nuc_mxtotal_offgrid = sum( wh07_nuc_mx(wh07_nuc_offgrid,:), 1 ) + wh07_nuc_mx_wind;
wh07_nuc_xntotal_offgrid = wh07_nuc_mxtotal_offgrid ./ sum( wh07_nuc_mxtotal_offgrid );

wh07_nuc_xnsolar                = zeros(size(xn_solar));
wh07_nuc_xnsolar(solar_in_wh07) = wh07_nuc_xntotal_ejecta(wh07_solar);
wh07_nuc_xn_pf                  = wh07_nuc_xnsolar ./ xn_solar;

wh07_nuc_ysolar                = zeros(size(xn_solar));
wh07_nuc_ysolar(solar_in_wh07) = wh07_nuc_ytotal_ejecta(wh07_solar);
wh07_nuc_y_pf                  = wh07_nuc_ysolar ./ y_solar;

wh07_nuc_xntotal_z_ejecta(1:max(wh07_zz)) = 0.0;
wh07_nuc_ytotal_z_ejecta(1:max(wh07_zz)) = 0.0;
for i = 1:max(wh07_zz)
    wh07_nuc_xntotal_z_ejecta(i) = sum( wh07_nuc_xntotal_ejecta(inuc_z_wh07(i)) );
    wh07_nuc_ytotal_z_ejecta(i)  = sum( wh07_nuc_ytotal_ejecta(inuc_z_wh07(i)) );
end

wh07_nuc_xntotal_a_ejecta(1:max(wh07_aa)) = 0.0;
wh07_nuc_ytotal_a_ejecta(1:max(wh07_aa)) = 0.0;
for i = 1:max(wh07_aa)
    wh07_nuc_xntotal_a_ejecta(i) = sum( wh07_nuc_xntotal_ejecta(inuc_a_wh07(i)) );
    wh07_nuc_ytotal_a_ejecta(i)  = sum( wh07_nuc_ytotal_ejecta(inuc_a_wh07(i)) );
end

wh07_nuc_xn_z_pf = zeros( size(xn_z_solar) );
wh07_nuc_y_z_pf = zeros( size(xn_z_solar) );
wh07_nuc_xn_a_pf = zeros( size(xn_a_solar) );
wh07_nuc_y_a_pf = zeros( size(xn_a_solar) );
wh07_nuc_xn_z_pf(1:min(max(wh07_zz),max(zz_solar))) = wh07_nuc_xntotal_z_ejecta(1:min(max(wh07_zz),max(zz_solar))) ./ xn_z_solar(1:min(max(wh07_zz),max(zz_solar))); 
wh07_nuc_y_z_pf(1:min(max(wh07_zz),max(zz_solar)))  = wh07_nuc_ytotal_z_ejecta(1:min(max(wh07_zz),max(zz_solar)))  ./ y_z_solar(1:min(max(wh07_zz),max(zz_solar)));
wh07_nuc_xn_a_pf(1:min(max(wh07_aa),max(aa_solar))) = wh07_nuc_xntotal_a_ejecta(1:min(max(wh07_aa),max(aa_solar))) ./ xn_a_solar(1:min(max(wh07_aa),max(aa_solar))); 
wh07_nuc_y_a_pf(1:min(max(wh07_aa),max(aa_solar)))  = wh07_nuc_ytotal_a_ejecta(1:min(max(wh07_aa),max(aa_solar)))  ./ y_a_solar(1:min(max(wh07_aa),max(aa_solar))); 

save( wh07_nuc_fname, ...
    'wh07_nuc_*' );

% 1 year data
wh07_fin_nspecies = wh07_nspecies;
[ ~, ~, wh07_fin_r, ~, wh07_fin_vr, wh07_fin_rho, ~, wh07_fin_temp, wh07_fin_press, wh07_fin_eint, wh07_fin_enpy, wh07_fin_vth, wh07_fin_abar, wh07_fin_ye, wh07_fin_xn, ~, ~, ~, ~ ] ...
    = import_wh07( wh07data_fin_fname, wh07_fin_nspecies, wh07_fin_start_line, wh07_fin_end_line );

wh07_fin_offgrid  = wh07_offgrid;
wh07_fin_m  = wh07_m;
wh07_fin_dm = wh07_dm;
wh07_fin_mx = wh07_fin_xn .* repmat( wh07_fin_dm, 1, wh07_fin_nspecies );

wh07_fin_mx_wind = import_wh07_wind( wh07data_finwind_fname, wh07_fin_nspecies );
wh07_fin_mxtotal_ejecta  = sum( wh07_fin_mx(2:end,:), 1 ) + wh07_fin_mx_wind;
wh07_fin_xntotal_ejecta  = wh07_fin_mxtotal_ejecta ./ sum( wh07_fin_mxtotal_ejecta );
wh07_fin_ytotal_ejecta   = wh07_fin_xntotal_ejecta ./ wh07_aa';
wh07_fin_mxtotal_offgrid = sum( wh07_fin_mx(wh07_fin_offgrid,:), 1 ) + wh07_fin_mx_wind;
wh07_fin_xntotal_offgrid = wh07_fin_mxtotal_offgrid ./ sum( wh07_fin_mxtotal_offgrid );

wh07_fin_xnsolar                = zeros(size(xn_solar));
wh07_fin_xnsolar(solar_in_wh07) = wh07_fin_xntotal_ejecta(wh07_solar);
wh07_fin_xn_pf                  = wh07_fin_xnsolar ./ xn_solar;

wh07_fin_ysolar                = zeros(size(xn_solar));
wh07_fin_ysolar(solar_in_wh07) = wh07_fin_ytotal_ejecta(wh07_solar);
wh07_fin_y_pf                  = wh07_fin_ysolar ./ y_solar;

wh07_fin_xntotal_z_ejecta(1:max(wh07_zz)) = 0.0;
wh07_fin_ytotal_z_ejecta(1:max(wh07_zz)) = 0.0;
for i = 1:max(wh07_zz)
    wh07_fin_xntotal_z_ejecta(i) = sum( wh07_fin_xntotal_ejecta(inuc_z_wh07(i)) );
    wh07_fin_ytotal_z_ejecta(i)  = sum( wh07_fin_ytotal_ejecta(inuc_z_wh07(i)) );
end

wh07_fin_xntotal_a_ejecta(1:max(wh07_aa)) = 0.0;
wh07_fin_ytotal_a_ejecta(1:max(wh07_aa)) = 0.0;
for i = 1:max(wh07_aa)
    wh07_fin_xntotal_a_ejecta(i) = sum( wh07_fin_xntotal_ejecta(inuc_a_wh07(i)) );
    wh07_fin_ytotal_a_ejecta(i)  = sum( wh07_fin_ytotal_ejecta(inuc_a_wh07(i)) );
end

wh07_fin_xn_z_pf = zeros( size(xn_z_solar) );
wh07_fin_y_z_pf = zeros( size(xn_z_solar) );
wh07_fin_xn_a_pf = zeros( size(xn_a_solar) );
wh07_fin_y_a_pf = zeros( size(xn_a_solar) );
wh07_fin_xn_z_pf(1:min(max(wh07_zz),max(zz_solar))) = wh07_fin_xntotal_z_ejecta(1:min(max(wh07_zz),max(zz_solar))) ./ xn_z_solar(1:min(max(wh07_zz),max(zz_solar))); 
wh07_fin_y_z_pf(1:min(max(wh07_zz),max(zz_solar)))  = wh07_fin_ytotal_z_ejecta(1:min(max(wh07_zz),max(zz_solar)))  ./ y_z_solar(1:min(max(wh07_zz),max(zz_solar)));
wh07_fin_xn_a_pf(1:min(max(wh07_aa),max(aa_solar))) = wh07_fin_xntotal_a_ejecta(1:min(max(wh07_aa),max(aa_solar))) ./ xn_a_solar(1:min(max(wh07_aa),max(aa_solar))); 
wh07_fin_y_a_pf(1:min(max(wh07_aa),max(aa_solar)))  = wh07_fin_ytotal_a_ejecta(1:min(max(wh07_aa),max(aa_solar)))  ./ y_a_solar(1:min(max(wh07_aa),max(aa_solar))); 

save( wh07_fin_fname, ...
    'wh07_fin_*' );

%% Initial radii of particle rows
prad_rows = import_prad_rows( prad_rows_fname );

r0           = reshape( repmat( prad_rows, 40, 1 ), 1, [] );
r0_unshocked = r0( 1, p_unshocked{end} );
m0_unshocked = interp1( wh07_presn_r.^3, wh07_presn_m, r0_unshocked.^3 );

save( pdata_fname, ...
    'prad_rows', '-append' );

%% Save MAT-files, scalars, and a few small arrays for easy access

temp_peak = max(temp);

ye_nse_8GK     = zeros(size(plist));
ye_nse_chimera = zeros(size(plist));
for p_id = plist
    ye_nse_8GK(p_id)     = ye(itime_nse_8GK(p_id),p_id);
    ye_nse_chimera(p_id) = ye(itime_nse_chimera(p_id),p_id);
end

save( model_fname, ...
    '*_fname','*_fname_base','model*','raw','p','pe','nu','pp','m','nse','shock','chim','grid','wh07','wh07_presn','wh07_nuc','wh07_fin','solar', ...
    'G','M_solar','c','dmnp','dunui','ergfoe','ergmev','h','kmev','ku','me','ncoef','rmu', ...
    'M_inner','M_tracer','M_masscut','M_pmasscut','r_cut','r_pcut','r_core','r_si','radius_inner','max_radius', ...
    'plist','p_grid','plist_all','p_extra','prad_rows','tracer_mass', ...
    'r0','r0_unshocked','m0_unshocked', ...
    'inuc_*','itime*','time*','i_*','iframe_*', ...
    'delta','rdelta', ...
    'nname_solar','aa_solar','zz_solar','xn_solar','y_solar','xn_a_solar','xn_z_solar','y_a_solar','y_z_solar', ...
    'wh07_aa','wh07_zz','wh07_solar','wh07_offgrid','wh07_nspecies','wh07_nname','wh07_m','wh07_dm','solar_wh07','solar_in_wh07','wh07_in_solar', ...
    'wh07_presn_offgrid','wh07_presn_nspecies', ...
    'temp_peak','ye_nse*' );

save( model_work_fname, ...
    '*_fname','*_fname_base','model*','raw','p','pe','nu','pp','m','nse','shock','chim','grid','wh07','wh07_presn','wh07_nuc','wh07_fin','solar', ...
    'G','M_solar','c','dmnp','dunui','ergfoe','ergmev','h','kmev','ku','me','ncoef','rmu', ...
    'M_inner','M_tracer','M_masscut','M_pmasscut','r_cut','r_pcut','r_core','r_si','radius_inner','max_radius', ...
    'plist','p_grid','plist_all','p_extra','prad_rows','tracer_mass', ...
    'r0','r0_unshocked','m0_unshocked', ...
    'inuc_*','itime*','time*','i_*','iframe_*', ...
    'delta','rdelta', ...
    'nname_solar','aa_solar','zz_solar','xn_solar','y_solar','xn_a_solar','xn_z_solar','y_a_solar','y_z_solar', ...
    'wh07_aa','wh07_zz','wh07_solar','wh07_offgrid','wh07_nspecies','wh07_nname','wh07_m','wh07_dm','solar_wh07','solar_in_wh07','wh07_in_solar', ...
    'wh07_presn_offgrid','wh07_presn_nspecies','wh07_presn_r','wh07_presn_m', ...
    'temp_peak','ye_nse*', ...
    'wh07_*_mx','wh07_*_xn', 'wh07_*_pf', ...
    'me_grid','rade_grid','rho_grid','ye_grid','chim_*','p_*', 'm_*' );