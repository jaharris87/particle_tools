%% fd_fit
% Fit neutrino spectrum to Fermi-Dirac distribution,
%
% $$f_{FD}(E_\nu;T_\nu,\eta) = \frac{1}{e^{\frac{E_\nu}{T_\nu} - \eta}},$$
%
% by matching the zeroth, first, and second energy moments.
%   i.e.
%
% $$\frac{\epsilon^{(0)} \epsilon^{(2)}}{(\epsilon^{(1)})^2} = \frac{F_4(\eta) F_2(\eta)}{(F_3(\eta))^2}
%                                                            = \frac{<E^2_\nu>}{<E_\nu>^2},$$
%
% where $F_n(\eta)$ is the complete Fermi-Dirac integral of order n and degeneracy $\eta$. (See Janka & Hillebrandt 1989 for details)
%
%% Function Definition
function [ eta, T_nu, residual ] = fd_fit( e0, e1, e2, varargin )
%% Input Arguments
%   Variable        Type        Dimension   Description
%   --------        -------     ---------   -------------------------------------------------
%   e0              Float       (M, N)      2D array of zeroth energy moment for N neutrino
%                                               types at each time interval
%
%   e1              Float       (M, N)      2D array of first energy moment for N neutrino
%                                               types at each time interval
%
%   e2              Float       (M, N)      2D array of second energy moment for N neutrino
%                                               types at each time interval
%
%% Output Arguments
%   Variable        Type        Dimension   Description
%   ---------       -------     ---------   -------------------------------------------------
%   eta             Float       (M, N)      Degeneracy parameter corresponding to fit
%   T_nu            Float       (M, N)      Neutrino temperature corresponding to fit
%   residual        Float       (M, N)      Absolute error in fit of the moments
%
%% Optional Input Arguments
%   Flag            Variable        Type        Dimension   Description
%   ---------       --------        -------     ---------   -------------------------------------------------
%   Tol             tol             Float        1          Tolerance criteria for iteration convergence
%     >=eps                                                     ( Default: 1.0E-05 )
%
%   MaxIter         itmax           Integer      1          Maximum number of iterations for convergence
%     >0                                                        ( Default: 20 )
%
%   PrintData       print           Logical      1          Flag for displaying report on iteration data
%                                                               ( Default: false )
%
%% Initialization

tic;

% Create an instance of the inputParser class.
p = inputParser;
p.KeepUnmatched = true;

% Define required inputs
p.addRequired('e0', ...
    @(x)validateattributes(x, {'numeric'}, ...
    {'2d', ...
    'real'}));
p.addRequired('e1', ...
    @(x)validateattributes(x, {'numeric'}, ...
    {'2d', ...
    'real'}));
p.addRequired('e2', ...
    @(x)validateattributes(x, {'numeric'}, ...
    {'2d', ...
    'real'}));

% Define optional inputs
p.addOptional('Tol', 1.0E-05, ...
    @(x)validateattributes(x, {'numeric'}, ...
    {'scalar', ...
    'real', ...
    '>=', eps}));
p.addOptional('MaxIter', 20, ...
    @(x)validateattributes(x, {'numeric'}, ...
    {'scalar', ...
    'integer', ...
    'positive'}));
p.addOptional('PrintData', false, ...
    @(x)validateattributes(x, {'logical'}, ...
    {'scalar'}));

% Parse, validate, and assign input arguments
p.parse( e0, e1, e2, varargin{:} );
e0    = p.Results.e0;
e1    = p.Results.e1;
e2    = p.Results.e2;
tol   = p.Results.Tol;
itmax = p.Results.MaxIter;
print = p.Results.PrintData;

% Get the number of neutrino types and timesteps from the dimensions of inputs
tsteps = size(e0,1);
nnu    = size(e0,2);

% Initialize output to non-degenerate values
eta   = zeros(tsteps,nnu);

%% Test that a solution exists
% If $\frac{<E^2_\nu>}{<E_\nu>^2}$ is larger than the corresponding moments from the F-D distriubtion
% at $\eta = 0$, then a solution doesn't exist; approximate the distribution as non-degenerate.

e2ee     = e0 .* e2 ./ e1.^2;
residual = e2ee - 486000 * zeta(3) * zeta(5) / ( 49 * pi^8 );
solution_exists = ( residual > 0 );

%% Fit the data

% Loop over each neutrino species
for n = 1:nnu
    
    % Loop over each timestep
    for t = 1:tsteps      
        
        if ~solution_exists(t,n)
            
            % Print status
            if print
                fprintf('\tFit does not exist for (t,n) = (%d,%d), assuming non-degenerate case\n',t,n);
            end
            
        else
            %% Initialize iteration variables
           
            % Store initial guess for each method to use
            etahold   = 0.1;
            converged = false;
            
            %% Secant Method
            
            % Initial values for first and second guesses of eta
            eta0 = etahold;
            eta1 = etahold+1;
            
            % Iterate
            for j = 1:itmax
                
                f2_0 = FD_int_approx(2,eta0);
                f3_0 = FD_int_approx(3,eta0);
                f4_0 = FD_int_approx(4,eta0);
                f_e2ee_0 = f2_0*f4_0/f3_0^2;
                dfeta0 = f_e2ee_0 - e2ee(t,n);
                
                f2_1 = FD_int_approx(2,eta1);
                f3_1 = FD_int_approx(3,eta1);
                f4_1 = FD_int_approx(4,eta1);
                f_e2ee_1 = f2_1*f4_1/f3_1^2;
                dfeta1 = f_e2ee_1 - e2ee(t,n);
                
                top = dfeta1*(eta1-eta0);
                bottom = dfeta1-dfeta0;
                
                % Update the solution
                etat = eta1 - top / bottom;
                
                % Check for convergence
                residualt = abs(etat-eta1);
                if residualt < tol
                    converged = true;
                    break
                else
                    % Update guess
                    eta0 = eta1;
                    eta1 = etat;
                end
            end
            
            %% Newton-Raphson
            % If Secant Method fails, try Newton-Raphson
            if ~converged
                
                % Restore initial guess
                eta0 = etahold;
                
                % Iterate
                for i = 1:itmax
                    
                    f1 = FD_int_approx(1,eta0);
                    f2 = FD_int_approx(2,eta0);
                    f3 = FD_int_approx(3,eta0);
                    f4 = FD_int_approx(4,eta0);
                    f_e2ee = f2*f4/f3^2;
                    
                    % Formulate f(eta) such that the solution corresponds to f(eta) = 0
                    f = f_e2ee - e2ee(t,n);
                    
                    % Calculate df(eta)/d(eta) ( d(F_n(eta))/d(eta) = F_{n-1}(eta) )
                    fprime = f_e2ee * ( f1/f2 + f3/f4 - 2*f2/f3 );
                    
                    % Update solution
                    etat = eta0 - f/fprime;
                    
                    % Check for convergence
                    residualt = abs(etat-eta0);
                    if residualt < tol
                        converged = true;
                        break
                    else
                        % Update guess
                        eta0 = etat;
                    end
                end
            end
            
            %% Bisection Method
            % If both Secant Method and Newton-Raphson fail, try Bisection Method
            if ~converged
                
                % Set bounds on the domain of eta
                etamax = 20.0;
                etamin = 0.0;
                
                f2_min = FD_int_approx(2,etamin);
                f3_min = FD_int_approx(3,etamin);
                f4_min = FD_int_approx(4,etamin);
                f_e2ee_min = f2_min*f4_min/f3_min^2;
                dfetamin = f_e2ee_min - e2ee(t,n);
                
                f2_max = FD_int_approx(2,etamax);
                f3_max = FD_int_approx(3,etamax);
                f4_max = FD_int_approx(4,etamax);
                f_e2ee_max = f2_max*f4_max/f3_max^2;
                dfetamax = f_e2ee_max - e2ee(t,n);
                
                % Check that there exists an eta between etamin and etamax such that f(eta) = 0
                if dfetamin*dfetamax < 0
                    
                    % Iterate
                    for k = 1:itmax
                        
                        % Bisect the interval as a guess to solution
                        eta = 0.5 * ( etamax + etamin );
                        
                        f2 = FD_int_approx(2,eta);
                        f3 = FD_int_approx(3,eta);
                        f4 = FD_int_approx(4,eta);
                        f_e2ee = f2*f4/f3^2;
                        dfeta = f_e2ee - e2ee(t,n);
                        
                        % Check for convergence
                        residualt = abs(dfeta);
                        if residualt < tol
                            converged = true;
                            break
                        else
                            % Update the bounds on the solution
                            if dfeta*dfetamin < 0
                                etamax = eta;
                                dfetamax = dfeta;
                            elseif dfeta*dfetamax < 0
                                etamin = eta;
                                dfetamin = dfeta;
                            end
                        end
                    end
                end
            end
            
            %% Finalize
            
            % If no solution was found, assume eta = 0 (non-degenerate)
            if ~converged
                if print
                    fprintf('\tSolution could not be found for (t,n) = (%d,%d); assuming non-degenerate case\n',t,n);
                end
            else
                eta(t,n) = etat;
                residual(t,n) = residualt;
            end
            
        end
    end    
end

%% Calculate the neutrino temperature
T_nu = e1 ./ e0 .* FD_int_approx(2,eta) ./ FD_int_approx(3,eta);

toc;

end

