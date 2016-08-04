%% get_iso
% Determine index for specified isotope
function [ inuc, inuc_aa, inuc_zz ] = get_iso( aa, zz, varargin )
%% Input Arguments
%   Variable        Type        Dimension   Units   Description
%   --------        -------     ---------   -----   ------------------------------------------------------------
%   aa              Integer      N                  Mass number (A) for each isotope
%     >0
% 
%   zz              Integer      N                  Atomic/proton number (Z) for each isotope
%     >=0
%
%% Output Arguments
%   Variable        Type        Dimension   Units   Description
%   ---------       -------     ---------   -----   ------------------------------------------------------------
%   inuc            Integer      M                  Indices of matching isotopes
%
%   inuc_aa         Integer      MA                 Indices of isotopes matching given A
%
%   inuc_zz         Integer      MZ                 Indices of isotopes matching given Z
%
%% Optional Input Arguments
%   Flag            Variable        Type        Dimension   Units   Description                                                     Default
%   ---------       --------        -------     ---------   -----   ------------------------------------------------------------    ---------------------------
%   A               aa_iso          Integer      NA                 Value(s) of A for which to determine index                       1:max(aa)
%    >0
%
%   Z               zz_iso          Integer      NZ                 Value(s) of Z for which to determine index                       0:max(zz)
%    >0
%
%   Name            name_iso        String       1                  Abbreviated name of element or particular isotope               ''
%                                                                       (e.g. 'Ni' or 'Ni56' or '56Ni')
%   
%% Initialization

load_constants;
el = build_element_symbol();
neutron   = {'n','nt','nt1','n1','n01','neut','neutron'};
proton    = {'p'};
deuterium = {'d'};
tritium   = {'t'};

% Create an instance of the inputParser class.
p = inputParser;
p.FunctionName = 'get_iso';
p.KeepUnmatched = true;

% Define required inputs
p.addRequired('aa', ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'vector', ...
                                         'integer', ...
                                         '>', 0}));
p.addRequired('zz', ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'vector', ...
                                         'integer', ...
                                         '>=', 0}));
% Define optional inputs
p.addOptional('A', [], ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'vector', ...
                                         'integer', ...
                                         '>', 0}));
p.addOptional('Z', [], ...
              @(x)validateattributes(x, {'numeric'}, ...
                                        {'vector', ...
                                         'integer', ...
                                         '>=', 0}));
p.addOptional('Name', '', ...
              @(x)validateattributes(x, {'char'}, ...
                                        {'row'}));

% Parse, validate, and assign input arguments.
p.parse( aa, zz, varargin{:} );
aa = p.Results.aa;
zz = p.Results.zz;

aa_iso   = p.Results.A;
zz_iso   = p.Results.Z;
name_iso = p.Results.Name;

% Extract Z and A values from element name
if ~isempty(name_iso)
    
    if strcmp( name_iso, neutron );
        name_iso = 'n1';
    elseif strcmp( name_iso, proton )
        name_iso = 'h1';
    elseif strcmp( name_iso, deuterium )
        name_iso = 'h2';
    elseif strcmp( name_iso, tritium )
        name_iso = 'h3';
    end
    
    [i,j] = regexp( name_iso, '[0-9]*' );
    if ~isempty(i) && ~isempty(j)
        aa_name = str2double(name_iso(i:j));
        aa_iso  = union( aa_iso, aa_name );
    end
    
    el_iso = regexprep( name_iso, '[0-9]*', '' );
    if strcmp( name_iso, 'n1' )
        zz_iso  = union( zz_iso, 0 );
    elseif ~isempty(el_iso)
        zz_name = find( strcmpi( cellstr(el), el_iso ) ) - 1;
        zz_iso  = union( zz_iso, zz_name );
    end
end

% Match all A if not specified
if isempty(aa_iso)
    aa_iso = 1:max(aa);
end

% Match all Z if not specified
if isempty(zz_iso)
    zz_iso = 0:max(zz);
end

% Make sure the isotope lists are in columns
aa     = reshape( aa, [], 1 );
zz     = reshape( zz, [], 1 );
aa_iso = reshape( aa_iso, [], 1 );
zz_iso = reshape( zz_iso, [], 1 );

% Find all isotopes which match given A
inuc_aa = find( ismember( aa, aa_iso ) );

% Find all isotopes which match given Z
inuc_zz = find( ismember( zz, zz_iso ) );

% Find all isotopes which match both A and Z
inuc = intersect( inuc_aa, inuc_zz );

% varargout{1} = inuc;
% varargout{2} = inuc_aa;
% varargout{3} = inuc_zz;

end