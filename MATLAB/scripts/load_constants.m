%% Physical constants
M_solar = 1.98892E+33;     % Solar mass [g]
G       = 6.67384E-08;     % Gravitational constant [cm^{3} g^{-1} s^{-2}]
c       = 2.99792458E+10;  % Speed of light [cm s^{-1}]
h       = 4.13567E-21;     % Planck's constant [MeV s]
hbarc   = 197.327E+00;     % hbar * c (MeV fm)
me      = 0.510998928E+00; % Electron mass [MeV c^{-2}]
kmev    = 8.61733E-11;     % Boltzmann's constant [MeV K^{-1}]
cm3fm3  = 1.0E-39;         % cm**3/fm**3
ergmev  = 1.6021773E-06;   % Ergs per MeV
ergfoe  = 1.0E-51;         % Conversion factor for ergs to foe
rmu     = 1.67378E-24;     % Mean baryon mass [g]
mb      = 938.919E+00;     % Baryon mass [MeV/c^2]
mp_ex   = 7.2889705E+00;   % proton mass excess [MeV]
mn_ex   = 8.0713171E+00;   % neutron mass excess [MeV]
amu     = 1.036427E-18;    % atomic mass unit [MeV]
amuc    = 931.49406121;    % atomic mass unit [MeV/c^2]
mevg    = 5.60958885E+26;  % [MeV] per [g c^{2}]
jev     = 1.602176565E-19; % Joules per eV

%% Derived constants
ku            = ergmev/rmu;       % ( # nucleons/gram )( erg/mev ) [g^{-1} erg MeV^{-1}]
kfm           = cm3fm3/rmu;       % # nucleons/gram )( cm3/fm3 ) [g^{-1} cm^{3} fm^{-3}]
hbar          = h / (2*pi);       % Planck's constant divided by 2*pi [MeV s];
dmnp_mev      = 1.29332E+00;      % Neutron-proton mass difference [(MeV)]
dmnp          = 1.29332E+00 * ku; % Neutron-proton mass difference [(MeV) * (g^{-1} erg MeV^{-1})]
pe_int_offset = 8.9E+00 * ku;     % Internal energy offset [(MeV) * (g^{-1} erg MeV^{-1})]
ncoef         = 4*pi / (h*c)^3;   % Neutrino distribution coefficient [ MeV^{-3} cm^{-3} ]

%% CHIMERA neutrino energy grid at infinity
unumn = 4;
unumx = 250;
nnugpmx = 20;
[unui,dunui] = nu_egrid(unumn, unumx, nnugpmx);
clear unumn unumx nnugpmx