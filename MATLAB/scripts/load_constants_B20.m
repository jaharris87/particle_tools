%% Model-specific constants
M_inner      = 2.83169E+33;     % Inner mass not captured by particles [g]
M_tracer     = 7.05069E+29;     % Individual tracer mass [g]
radius_inner = 1.38533E+08;     % Inner radius corresponding to inner mass [cm]
max_radius   = 1.99553708E+09;  % Radius of outer-most zone in initial model [cm]
time_bounce  = 4.60602308E-01;  % Time of core bounce [s]
time_final   = 2.15278;         % Final time shared by all tracers [s]
time_bseries = 1.843602308;     % End time used in B-series paper [s]

M_masscut    = 3.6284570993849653E+33; % Mass-cut from 1D WH07 model [g]
r_core       = 1.5907993162120250E+08; % Outer radius of inner iron-core [cm]
r_si         = 2.6905311727127808E+08; % Outer radius of silicon shell [cm]

% Chimera frame data
ny_frame    = 256;
start_frame = 1;
end_frame   = 8562;

%% Physical constants
load_constants;