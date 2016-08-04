function y = FD_int_approx( j, eta )

% Analytic approximations for Fermi-Dirac integrals of order j > -1/2
% Date: September 29, 2008
% Author: Raseong Kim (Purdue University)
%
% Inputs
% eta: eta_F
% j: FD integral order
%
% Outputs
% y: value of FD integral (the "script F" defined by Blakemore (1982))
%
% For more information in Fermi-Dirac integrals, see:
% "Notes on Fermi-Dirac Integrals (3rd Edition)" by Raseong Kim and Mark
% Lundstrom at http://nanohub.org/resources/5475
% 
% References
% [1]D. Bednarczyk and J. Bednarczyk, Phys. Lett. A, 64, 409 (1978)
% [2]J. S. Blakemore, Solid-St. Electron, 25, 1067 (1982)
% [3]X. Aymerich-Humet, F. Serra-Mestres, and J. Millan, Solid-St. Electron, 24, 981 (1981)
% [4]X. Aymerich-Humet, F. Serra-Mestres, and J. Millan, J. Appl. Phys., 54, 2850 (1983)

if j < -1/2
    error( 'The order should be equal to or larger than -1/2.')
else
    x = eta;
    switch j
        case 0
            y = log( 1 + exp( x ) );        % analytic expression
            
        case 1/2
            % Model proposed in [1]
            % Expressions from eqs. (22)-(24) of [2]
            mu = x .^ 4 + 50 + 33.6 * x .* ( 1 - 0.68 * exp( -0.17 * ( x + 1 ) .^ 2 ) );
            xi = 3 * sqrt( pi ) ./ ( 4 * mu .^ ( 3 / 8 ) );
            y = ( exp( - x ) + xi ) .^ -1;
        
        case 3/2
            % Model proposed in [3]
            % Expressions from eq. (5) of [3]
            % The integral is divided by gamma( j + 1 ) to make it consistent with [1] and [2].
            a = 14.9;
            b = 2.64;
            c = 9 / 4;
            y = ( ( j + 1 ) * 2 ^ ( j + 1 ) ./ ( b + x + ( abs( x - b ) .^ c + a ) .^ ( 1 / c ) ) .^ ( j + 1 ) ...
            + exp( -x ) ./ gamma( j + 1 ) ) .^ -1 ./ gamma( j + 1 );
            
        otherwise
            % Model proposed in [4]
            % Expressions from eqs. (6)-(7) of [4]
            % The integral is divided by gamma( j + 1 ) to make it consistent with [1] and [2].
            a = ( 1 + 15 / 4 * ( j + 1 ) + 1 / 40 * ( j + 1 ) ^ 2 ) ^ ( 1 / 2 );
            b = 1.8 + 0.61 * j;
            c = 2 + ( 2 - sqrt( 2 ) ) * 2 ^ ( - j );
            y = ( ( j + 1 ) * 2 ^ ( j + 1 ) ./ ( b + x + ( abs( x - b ) .^ c + a ^ c ) .^ ( 1 / c ) ) .^ ( j + 1 ) ...
                + exp( -x ) ./ gamma( j + 1 ) ) .^ -1;
    end
end
