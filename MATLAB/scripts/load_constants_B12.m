%% Model-specific constants
M_inner      = 2.39284E+33;     % Inner mass not captured by particles [g]
M_tracer     = 3.71525E+29;     % Individual tracer mass [g]
radius_inner = 8.90841E+07;     % Inner radius corresponding to inner mass [cm]
max_radius   = 2.90390188E+09;  % Radius of outer-most zone in initial model [cm]
time_bounce  = 2.63176547E-01;  % Time of core bounce [s]
time_final   = 1.67357;         % Final time shared by all tracers [s]
time_bseries = 1.599176547;     % End time used in B-series paper [s]

M_masscut    = 3.0434759993825662E+33; % Mass-cut from 1D WH07 model [g]
r_core       = 1.1217214652561970E+08; % Outer radius of inner iron-core [cm]
r_si         = 2.8132029403125674E+08; % Outer radius of silicon shell [cm]

% Chimera frame data
ny_frame    = 256;
start_frame = 1;
end_frame   = 7122;

%% Physical constants
load_constants;