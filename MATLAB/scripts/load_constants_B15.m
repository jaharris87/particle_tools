%% Model-specific constants
M_inner      = 2.62787E+33;     % Inner mass not captured by particles [g]
M_tracer     = 5.69632E+29;     % Individual tracer mass [g]
radius_inner = 1.13129E+08;     % Inner radius corresponding to inner mass [cm]
max_radius   = 1.96231147E+09;  % Radius of outer-most zone in initial model [cm]
time_bounce  = 3.3247634E-01;   % Time of core bounce [s]
time_final   = 2.00218;         % Final time shared by all tracers [s]
time_bseries = 1.53247634;      % End time used in B-series paper [s]

M_masscut    = 3.6521711997476101E+33; % Mass-cut from 1D WH07 model [g]
r_core       = 1.3999434804229215E+08; % Outer radius of inner iron-core [cm]
r_si         = 3.9377265840726244E+08; % Outer radius of silicon shell [cm]

% Chimera frame data
ny_frame    = 256;
start_frame = 1;
end_frame   = 8415;

%% Physical constants
load_constants;