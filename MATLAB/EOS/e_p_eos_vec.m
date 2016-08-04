function [ pe_out, ee_out, se_out, ue_out, yeplus_out, rel_out, rel1_out, rel2_out ] = e_p_eos_vec( brydns_in, t_mev_in, ye_in )

brydns = reshape(brydns_in,1,[]);
t_mev = reshape(t_mev_in,1,[]);
ye = reshape(ye_in,1,[]);

nx = length(t_mev);
pe(1:nx) = 0.0;
ee(1:nx) = 0.0;
se(1:nx) = 0.0;
ue(1:nx) = 0.0;
yeplus(1:nx) = 0.0;
rel(1:nx) = 0.0;
rel1(1:nx) = 0.0;
rel2(1:nx) = 0.0;

zero = 0.0;
third = 1.0/3.0;

hbarc   = 197.33E+00;      % hbar * c (GeV fm)
me      = 0.510998928E+00; % Electron mass [MeV c^{-2}]

nlag = 48;
ncnvge = 30;

eps1 = 1.0e-4;
eps2 = 1.0e-6;
beta_max = 2.0;
relmin = 1.2;
etabet = 2.0;
tthird = 2.0/3.0;

[ xa, wt0 ] = glaquad( nlag );

wta = fexp( xa ) .* wt0;

pi2 = pi^2;
pi4 = pi2^2;
const = hbarc^3 * pi2;
cnstnr = hbarc^3 * 3.0 * pi2 / sqrt( 2.0 * me^3 );
coef_e_fermi_nr = ( hbarc^2 / ( 2.0 * me ) ) * ( 3.0 * pi2 )^tthird;

rmuec = 0.5 * ( 2.0 * pi * hbarc^2 / me )^1.5;

e_fermi_nr = coef_e_fermi_nr * ( brydns .* ye ).^tthird;
non_rel(1:nx) = false;
non_rel( e_fermi_nr < 0.01 * me & t_mev < 0.01 * me ) = true;

ne_x_coef(1:nx)     = zero;
pfc(1:nx)           = zero;
beta(1:nx)          = zero;
etae(1:nx)          = zero;
ped(1:nx)           = zero;
tolp(1:nx)          = zero;
Failure(1:nx)       = false;
iterate(1:nx)       = true;
wwchk(1:nx)         = zero;
deriv(1:nx)         = zero;
fn_a(1:nlag,1:nx)   = zero;
fp_a(1:nlag,1:nx)   = zero;
fact(1:nlag,1:nx)   = zero;
tol(1:nx)           = Inf;
deta(1:nx)          = zero;
detap(1:nx)         = zero;
high_t(1:nx)        = false;
use_high_t(1:nx)    = false;
rel_gl(1:nx)        = false;

ne_x_coef(~non_rel) = const  * brydns(~non_rel)                ./ t_mev(~non_rel).^3.0;
ne_x_coef( non_rel) = cnstnr * brydns( non_rel) .* ye(non_rel) ./ t_mev( non_rel).^1.5;

pfc                 = hbarc * ( 3.0 * pi2 * brydns .* ye ).^third;

beta                = me ./ t_mev;
beta2               = beta .* beta;

etae                = cube( 1.5 * ( ne_x_coef .* ye ), max( pi2 * third - 0.5 * beta2, zero ) );
etae2               = etae .* etae;

high_t              = ( beta <= beta_max ) & ( beta <= tthird | etae > etabet * beta );
ped(high_t)         = ( t_mev(high_t) ./ ne_x_coef(high_t) ) * third .* ...
                        ( g3( etae2(high_t) ) - 1.5 * beta2(high_t) .* g1( etae2(high_t) ) );
pe(high_t)          = brydns(high_t) .* ped(high_t);
se(high_t)          = ( 4.0 ./ t_mev(high_t) ) .* ped(high_t) - ye(high_t) .* etae(high_t) + ...
                        beta2(high_t) .* g1( etae2(high_t) ) ./ ne_x_coef(high_t);
ee(high_t)          = t_mev(high_t) .* ( se(high_t) + ye(high_t) .* etae(high_t) ) - ped(high_t);
yeplus(high_t)      = 2.0 * fexp(-etae(high_t)) .* ( 1.0 - fexp(-etae(high_t)) / 8.0 ) ./ ne_x_coef(high_t);

rel1(high_t)        = 3.0d0 * ped(high_t) ./ ( ee(high_t) - me * ( ye(high_t) + 2.0 * yeplus(high_t) ) );
use_high_t          = ( rel1 <= relmin & rel1 ~= zero );

ue(use_high_t)      = etae(use_high_t) .* t_mev(use_high_t);
rel(use_high_t)     = rel1(use_high_t);

xa                  = repmat(xa',1,nx);
wta                 = repmat(wta',1,nx);

for j = 1:nlag
    fact(j,:)       = wta(j,:) .* sqrt( xa(j,:) .* ( xa(j,:) + 2.0 * beta ) ) .* ( xa(j,:) + beta );
end

for it = 1:ncnvge
    for j = 1:nlag
        fn_a(j,iterate)    = 1.0 ./ ( 1.0 + fexp( xa(j,iterate) + beta(iterate) - etae(iterate) ) );
        fp_a(j,iterate)    = 1.0 ./ ( 1.0 + fexp( xa(j,iterate) + beta(iterate) + etae(iterate) ) );
        wwchk(iterate)     = wwchk(iterate) + fact(j,iterate) .* ( fn_a(j,iterate) - fp_a(j,iterate) );
        deriv(iterate)     = deriv(iterate) + fact(j,iterate) .* ( fn_a(j,iterate) .* ...
                               ( 1.0 - fn_a(j,iterate) ) + fp_a(j,iterate) .* ( 1.0 - fp_a(j,iterate) ) );
    end
    Failure(wwchk <= zero) = true;
    iterate(Failure)       = false;
    
    tol(iterate)                = log( wwchk(iterate) ./ ( ne_x_coef(iterate) .* ye(iterate) ) );
    iterate(abs(tol) <= eps1)   = false;
    
    deriv(iterate)         = deriv(iterate) ./ wwchk(iterate);
    
    deta(tol .* tolp >= zero)   = -tol(tol .* tolp >= zero) ./ deriv(tol .* tolp >= zero);
    deta(tol .* tolp  < zero)   = -tol(tol .* tolp  < zero) .* detap(tol .* tolp  < zero) ./ ...
                                            ( tolp(tol .* tolp  < zero) - tol(tol .* tolp  < zero) );
                                        
    iterate(abs(deta) <= eps2 * abs(etae)) = false;
    
    deta(etae+deta <= zero)     = deta(etae+deta <= zero) / 10.0;
    
    etae(iterate)           = etae(iterate) + deta(iterate);
    
    etae(etae <= zero)      = abs(deta(etae <= zero));
    tolp(iterate)           = tol(iterate);
    detap(iterate)          = deta(iterate);
    
    if all(~iterate)
        break
    end
end

% Failure(iterate)            = true;

rel_gl                      = ~use_high_t & ~Failure & ( etae < 35.0 | pfc/me < 10.0 );

for j = 1:nlag
    fact(j,:)       = wta(j,:) .* sqrt( xa(j,:) .* ( xa(j,:) + 2.0 * beta ) );
end

beta_tmp                    = repmat(beta,nlag,1);

tmp                         = cumsum( fact .* ( fn_a + fp_a ) .* xa .* ( xa + 2.0 * beta_tmp ), 2 );
pe(rel_gl)                  = tmp(rel_gl) .* ( t_mev(rel_gl) ./ ne_x_coef(rel_gl) ) .* brydns(rel_gl) * third;

tmp                         = cumsum( fact .* ( fn_a + fp_a ) .* ( xa + beta_tmp).^2, 2 );
ee(rel_gl)                  = tmp(rel_gl) .* ( t_mev(rel_gl) ./ ne_x_coef(rel_gl) );
                                  
tmp                         = cumsum( fact .* fp_a .* ( xa + beta_tmp ), 2 );
yeplus(rel_gl)              = tmp(rel_gl) ./ ne_x_coef(rel_gl);

ue(rel_gl)                  = etae(rel_gl) .* t_mev(rel_gl);
se(rel_gl)                  = ( ee(rel_gl) + pe(rel_gl) ./ brydns(rel_gl) ) ./ t_mev(rel_gl) - ye(rel_gl) .* etae(rel_gl);

rel2(rel_gl)                = 3.0 * pe(rel_gl) ./ brydns(rel_gl) ./ ( ee(rel_gl) - me * ( ye(rel_gl) + 2.0 * yeplus(rel_gl) ) );
rel(rel_gl)                 = rel2(rel_gl);

pe_out = reshape(pe,size(t_mev_in));
ee_out = reshape(ee,size(t_mev_in));
se_out = reshape(se,size(t_mev_in));
ue_out = reshape(ue,size(t_mev_in));
yeplus_out = reshape(yeplus,size(t_mev_in));
rel_out = reshape(rel,size(t_mev_in));
rel1_out = reshape(rel1,size(t_mev_in));
rel2_out = reshape(rel2,size(t_mev_in));

end