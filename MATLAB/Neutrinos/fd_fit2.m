function [ beta ] = fd_fit2( e0, e1, e2, beta0, beta1, tol, itmax )
%UNTITLED5 Summary of this function goes here
%   Dbetailed explanation goes here

betahold = beta0;
e2ee = e2 * e0 / e1^2;
i = itmax;
j = itmax;
k = itmax;

%% Newton-Raphson
for i = 1:itmax
    f = dfbeta(e0,e1,e2,beta0);
    fprime = dfbeta_dbeta(beta0);
    beta = beta0 - f/fprime;
    
    if abs(beta-beta0) < tol
        break
    else
        beta0 = beta;
    end
end

%% Secant Method
if i == itmax
    beta0 = betahold;
    for j = 1:itmax
        top = dfbeta(e0,e1,e2,beta1)*(beta1-beta0);
        bottom = dfbeta(e0,e1,e2,beta1)-dfbeta(e0,e1,e2,beta0);
        beta = beta1 - top / bottom;
        if abs(beta-beta1) < tol
            break
        else
            beta0 = beta1;
            beta1 = beta;
        end
    end
end     

%% Bisection Method
if j == itmax
    beta0 = betahold;
    if dfbeta(e0,e1,e2,beta0)*dfbeta(e0,e1,e2,beta1) < 0
        betamax = max(beta0,beta1);
        betamin = min(beta0,beta1);
        for k = 1:itmax
            beta = 0.5 * ( betamax + betamin );
            f = dfbeta(e0,e1,e2,beta);
            if f < tol
                break
            else
                if f*dfbeta(e0,e1,e2,betamin) < 0
                    betamax = beta;
                elseif f*dfbeta(e0,e1,e2,betamax) < 0
                    betamin = beta;
                end
            end
        end
    end
end

display([' beta = ',beta,'  abs. error: ',abs(g_e2ee(beta)-e2ee)]);
display([' NR iters.: ',i,'  Secant iters.: ',j,'  Bisection iters.: ',k]);
end

