%% Model-specific constants
M_inner      = 2.94404E+33;     % Inner mass not captured by particles [g]
M_tracer     = 6.93238E+29;     % Individual tracer mass [g]
radius_inner = 1.47233E+08;     % Inner radius corresponding to inner mass [cm]
max_radius   = 2.20091438E+09;  % Radius of outer-most zone in initial model [cm]
time_bounce  = 4.68425713E-01;  % Time of core bounce [s]
time_final   = 2.16725;         % Final time shared by all tracers [s]
time_bseries = 1.867425713;     % End time used in B-series paper [s]

M_masscut    = 3.7808106869243770E+33; % Mass-cut from 1D WH07 model [g]
r_core       = 1.6858016534685987E+08; % Outer radius of inner iron-core [cm]
r_si         = 2.8819201239934283E+08; % Outer radius of silicon shell [cm]

% Chimera frame data
ny_frame    = 256;
start_frame = 1;
end_frame   = 8587;


%% Physical constants
load_constants;