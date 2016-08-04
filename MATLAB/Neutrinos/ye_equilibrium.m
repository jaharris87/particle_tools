function [ye_eq0,ye_eq1] = ye_equilibrium(lum_nu,tmev_nu,lum_nubar,tmev_nubar,radius,tmev,etae)

me = 0.510998928E+00;
c = 2.99792458E+10;
ergmev = 1.6021773E-06;
ergfoe = 1e-51;

r7 = radius .* 1e-7;

c1 = 1 + 0.6283./tmev_nu + 0.1293./tmev_nu.^2;
c2 = 1 + 1.158./tmev_nubar + 0.6./tmev_nubar.^2 + 0.1778./tmev_nubar.^3;

lambda_nue = 0.1945 .* (lum_nu.*ergfoe) .* tmev_nu ./ r7.^2 .* c1;
lambda_nuebar = 0.2 .* exp(-1.804./tmev_nu) .* (lum_nubar*ergfoe) .* tmev_nubar ./ r7.^2 .* c2;

j1 = 1 + 1.158./tmev_nubar - 0.6283./tmev_nu;
lambda_ratio1 = 1.029 * (lum_nubar./lum_nu) .* (tmev_nubar./tmev_nu) .* exp( -1.804./tmev_nubar) .* j1;
% lambda_ratio1 = lambda_nuebar ./ lambda_nue;

ye_eq0 = 1 ./ (1 + lambda_ratio1);

c3 = 1 + 0.646./tmev + 0.128./tmev.^2;
c4 = 1 + 1.16./tmev + 0.601./tmev.^2 + 0.178./tmev.^3 + 0.035./tmev.^4;

lambda_e = 1.578e-2 * (tmev/(me*c^2)).^5 .* exp( (-1.293 + etae)./tmev ) .* c3;
lambda_ep = 1.578e-2 * (tmev/(me*c^2)).^5 .* exp( (-0.511 - etae)./tmev ) .* c4;

lambda_ratio2 = exp( (-0.782 + 2*etae)./tmev ) .* (c3./c4);
% lambda_ratio2 = lambda_e ./ lambda_ep;

lambda_ratio3 = 2.34 * tmev.^5 .* r7.^2 ./ lum_nu ./ tmev_nu .* exp( (-me - etae)./tmev ) .* c4 ./ c1;
% lambda_ratio3 = lambda_ep ./ lambda_nue;

ye_eq1 = ye_eq0.*(1+lambda_ratio3.*(1-(1+lambda_ratio2)./(1+lambda_ratio1)));

% ye_eq2 = ye_eq0 ...
%        + ye_eq0*lambda_ratio3*(1-(1+lambda_ratio2)/(1+lambda_ratio1)) ...
%        + x_alpha*(0.5-ye_eq0) ...
%        + 0.5*ye_eq0*x_alpha*lambda_ratio3 ...
%          * (lambda_ratio2-1-(lambda_ratio1-1)*((lambda_ratio2+1)/(lambda_ratio1+1)));
end