function [ F ] = F_eta( n, eta )
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

tol = 1e-5;
itmax = 100;
F = 0;

if n == 0
    
    F = log(1 + exp(eta));
    
elseif n == -1
    
    F = 1 / (1 + exp(-eta));
    
elseif n <= -2
    
    coef = zeros(-(n+2),-(n+2));
    coef(1:-(n+2),1) = 1;
    for k = 2:-(n+2)
        for i = k:-(n+2)
            coef(k,i) = i*coef(k-1,i) - (k-i+1)*coef(k-1,i-1);
        end
    end
    
    P = exp(eta);
    for i = 2:-(n+2)
        P = P + coef(-(n+2),i)*exp((i-1)*eta);
    end
    
    F = exp(eta) * P / (1 + exp(eta) )^(-n);

elseif eta == 0
    
    F = (1-2^(-n))*zeta(n+1);
    
elseif eta < 0
    
    F_k = exp(eta);
    k = 2;
    while abs(F_k) > tol && k-1 <= itmax
        F = F + F_k;
        F_k = (-1)^(k-1) * exp(k*eta) / k^(n+1);
        k = k + 1;
    end
    
elseif eta > 0
    
    F1 = 0;
    F1_j = 0.5 / gamma(n+2);
    j = 1;
    while abs(F1_j) > tol && j <= itmax
        F1 = F1 + F1_j;        
        F1_j = (1 - 2^(1-2*j)) * zeta(2*j) / gamma(n+2-2*j) / eta^(2*j);
        j = j + 1;
    end
    
    F2 = 0;
    F2_k = exp(-eta);
    k = 2;
    while abs(F2_k) > tol && k-1 <= itmax
        F2 = F2 + F2_k;
        F2_k = (-1)^(k-1) * exp(-k*eta) / k^(n+1);
        k = k + 1;
    end
    
    F = 2*eta^(n+1)*F1 + cos(pi*n)*F2;

end

F = F * gamma(n+1);

end
