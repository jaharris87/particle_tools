function [ pe_out, ee_out, se_out, ue_out, yeplus_out, rel_out, rel1_out, rel2_out, approx_out ] = e_p_eos( brydns_in, t_mev_in, ye_in )

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
approx_out(1:nx) = 0.0;

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

for i = 1:nx
    
    if ~non_rel(i)
        
        ne_x_coef     = brydns(i) .* const./t_mev(i).^3;
        pfc           = hbarc * ( 3.0 * pi2 * brydns(i) .* ye(i) ).^third;
        beta          = me ./ t_mev(i);
        beta2         = beta * beta;
        etae          = zero;

        approx = ( beta <= beta_max );
%         approx = false;
%         approx = true;
        
        if approx
            etae        = cube( 1.5 * ( ne_x_coef .* ye(i) ), max( pi2 * third - 0.5 * beta2, zero ) );
            approx      = ( beta <= tthird ) || ( etae > etabet * beta );
%             approx      = true;
            
            if approx
                ped       = ( t_mev(i)/ne_x_coef ) * third * ( g3( etae * etae ) - 1.5 * beta2 * g1( etae * etae ) );
                pe(i)     = brydns(i) * ped;
                se(i)     = ( 4.0 ./ t_mev(i) ) * ped - ye(i) * etae + beta2 * g1( etae * etae ) / ne_x_coef;
                ee(i)     = t_mev(i) * ( se(i) + ye(i) * etae ) - ped;
                yeplus(i) = 2.0 * fexp(-etae) * ( 1.0 - fexp(-etae) / 8.0 ) / ne_x_coef;
                
                rel1(i)    = 3.0d0 * ped ./ ( ee(i) - me * ( ye(i) + 2.0 * yeplus(i) ) );
                approx = ( rel1(i) <= relmin );
%                 approx = true;
            end
        end
        
        if approx
            ue(i) = etae * t_mev(i);
            rel(i) = rel1(i);
            approx_out(i) = 1;
            continue
        end
        
        etae = cube( 1.5 * ( ne_x_coef * ye(i) ) , max( pi2 * third - 0.5 * beta2, zero ) );
        
        tolp = zero;
        Failure = false;
        
        for it = 1:ncnvge
            wwchk = zero;
            deriv = zero;
            
            fn_a(1:nlag) = zero;
            fp_a(1:nlag) = zero;
            
            for j = 1:nlag
                fact      = wta(j) .* sqrt( xa(j) .* ( xa(j) + 2.0 * beta ) ) .* ( xa(j) + beta );
                fn_a(j)   = 1.0 ./ ( 1.0 + fexp( xa(j) + beta - etae ) );
                fp_a(j)   = 1.0 ./ ( 1.0 + fexp( xa(j) + beta + etae ) );
                wwchk     = wwchk + fact .* ( fn_a(j) - fp_a(j) );
                deriv     = deriv + fact .* ( fn_a(j) .* ( 1.0 - fn_a(j) ) + fp_a(j) .* ( 1.0 - fp_a(j) ) );
            end
            
            if wwchk <= zero
                Failure = true;
                break
            end
            
            tol = log( wwchk / ( ne_x_coef * ye(i) ) );
            if abs(tol) <= eps1
                break
            end
            deriv = deriv/wwchk;
            
            if tol * tolp >= zero
                deta = -tol/deriv;
            else
                deta = detap * tol / ( tolp - tol );
            end
            
            if abs(deta) <= eps2 * abs(etae)
                break
            end
            
            if etae + deta <= zero
                deta = deta / 10.0;
            end
            etae = etae + deta;
            
            if etae <= zero
                etae = abs(deta);
            end
            tolp = tol;
            detap = deta;
            
            if it == ncnvge
                Failure = true;
            end
        end
        
        if ~Failure && ( etae < 35.0 || pfc/me < 10.0 )
            pe(i) = zero;
            ee(i) = zero;
            se(i) = zero;
            yeplus(i) = zero;
            
            for j = 1:nlag
                fact      = wta(j) .* sqrt( xa(j) .* ( xa(j) + 2.0 * beta ) );
                pe(i)     = pe(i) + fact .* ( fn_a(j) + fp_a(j) ) .* xa(j) .* ( xa(j) + 2.0 * beta );
                ee(i)     = ee(i) + fact .* ( fn_a(j) + fp_a(j) ) .* ( xa(j) + beta ).^2;
                yeplus(i) = yeplus(i) + fact .* fp_a(j) .* ( xa(j) + beta );
            end
            
            yeplus(i) = yeplus(i)/ne_x_coef;
            ue(i) = etae * t_mev(i);
            pe(i) = ( t_mev(i) / ne_x_coef ) .* brydns(i) * third .* pe(i);
            ee(i) = ( t_mev(i) / ne_x_coef ) .* ee(i);
            se(i) = ( ee(i) + pe(i) ./ brydns(i) ) ./ t_mev(i) - ye(i) * etae;
            rel2(i) = 3.0 * pe(i) ./ brydns(i) ./ ( ee(i) - me * ( ye(i) + 2.0 * yeplus(i) ) );
            rel(i) = rel2(i);
            approx_out(i) = 2;
            continue
        end
        
        rb2           = pi2/beta2;
        xans          = ( pfc/me )^3;
        xx            = pfc/me;
        dxxp          = zero;
        
        Failure = false;
        for l = 1:50
            tol = ffn(xx,beta) - xans / xans;
            if abs(tol) <= 1.0e-7
                break
            end
            dxx = -tol * xans / dfndx(xx,beta);
            if dxx*dxxp < zero
                dxx = dxx / 2.0;
            end
            dxxp = dxx;
            xx = xx + dxx;
            if l == 50
                Failure = true;
            end
        end
        
        if ~Failure
            pe(i)       = me^4 / ( 24.0 * pi2 * hbarc^3 ) * ( f(xx) + ...
                          4.0 * rb2 * xx * sqrt( 1.0 + xx^2 ) + ...
                          7.0/15.0 * rb2^2 * sqrt( 1.0 + xx * xx ) * ( 2.0 * xx^2 - 1.0 )/xx^3 );
            ek          = me^4 / ( 24.0 * pi2 * hbarc^3 ) * ( gg(xx) + 4.0 * rb2 * ( sqrt( 1.0 + xx^2 ) * ...
                          ( 3.0 * xx^2 + 1.0 )/xx - ( 2.0 * xx^2 + 1.0 )/xx ) );
            ee(i)       = ek/brydns(i) + ye(i) * me;
            ue(i)       = me * sqrt( 1.0 + xx^2 );
            se(i)       = ( ee(i) + pe(i)/brydns(i) - ue(i) * ye(i) ) / t_mev(i);
            rel(i)      = 3.0 * pe(i)/ek;
            yeplus(i)   = zero;
            approx_out(i) = 3;
            
            continue
        end

        disp('FAILURE');
        continue
    else
        tolp          = zero;
        etad          = e_fermi_nr(i) ./ t_mev(i);
        
        Sommerfeld    = false;
        if e_fermi_nr(i) / t_mev(i) > 35
            Sommerfeld = true;
        end
        
        if ~Sommerfeld
            
            rmuend      = t_mev(i) .* log( brydns(i) .* ye(i) .* rmuec ./ t_mev(i).^1.5 );
            etand       = rmuend ./ t_mev(i);
            rintrp      = etad / ( 1.0 + etad );
            etae        = etand + rintrp * ( etad - etand );
            ne_x_coef   = cnstnr * brydns(i) .* ye(i) ./ t_mev(i).^1.5;
            
            Failure = false;
            for it = 1:ncnvge
                wwchk = zero;
                deriv = zero;
                
                fn_a(1:nlag) = zero;
                
                for j = 1:nlag
                    fact      = wta(j) .* sqrt( xa(j) );
                    fn_a(j)   = 1.0 ./ ( 1.0 + fexp( xa(j) - etae ) );
                    wwchk     = wwchk + fact .* fn_a(j);
                    deriv     = deriv + fact .* ( fn_a(j) .* ( 1.0 - fn_a(j) ) );
                end
                
                tol = log( wwchk / ne_x_coef );
                if abs(tol) <= eps1
                    break
                end
                deriv = deriv/wwchk;
                
                if tol * tolp >= zero
                    deta = -tol/deriv;
                else
                    deta = detap * tol / ( tolp - tol );
                end
                
                if abs(deta) <= eps2 * abs(etae)
                    break
                end
                
                etae = etae + deta;
                tolp = tol;
                detap = deta;
                
                if it == ncnvge
                    Failure = true;
                end
            end
            
            if ~Failure
                pe(i) = zero;
                ek    = zero;
                se(i) = zero;
                for j = 1:nlag
                    fact = wta(j) * xa(j)^1.5;
                    ek = ek + fact * fn_a(j);
                end
                ue(i) = etae * t_mev(i) + me;
                ek    = ek * t_mev(i) * ye(i) / ne_x_coef;
                ee(i) = ek + ye(i) * me;
                pe(i) = 2.0 * ek * brydns(i) / 3.0;
                se(i) = ( ee(i) + pe(i) / brydns(i) ) ./ t_mev(i) - ye(i) * ue(i) ./ t_mev(i);
                yeplus(i) = zero;
                rel(i) = zero;
                approx_out(i) = 4;
                continue
            end
        end
        
        if Sommerfeld || Failure
            
            yans        = 3.0 * pi2 * brydns(i) * ye(i) * ( hbarc^2/( 2.0 * me * t_mev(i) ) ).^1.5;
            y           = yans.^(tthird);
            dyp         = zero;
            
            Failure = false;
            for l = 1:50
                tol = ( ffnr(y) - yans ) / yans;
                if abs(tol) <= 1.0e-7
                    break
                end
                dy = -tol * yans / dfnrdy(y);
                if dy * dyp < zero
                    dy = dy / 2.0;
                end
                dyp = dy;
                y = y + dy;
                if l == 50
                    Failure = true;
                end
            end
            
            if ~Failure
                pe(i) = 2.0 * t_mev(i) * pnr(y)/( 15.0 * pi2 ) * ( 2.0 * me * t_mev(i)/( hbarc^2 ) )^1.5;
                ek    = 1.5 * pe(i);
                ee(i) = ek ./ brydns(i) + ye(i) * me;
                ue(i) = y * t_mev(i) + me;
                se(i) = ( ee(i) + pe(i)./brydns(i) - ue(i) * ye(i) ) ./ t_mev(i);
                yeplus(i) = zero;
                rel(i) = zero;
                approx_out(i) = 5;
            end
        end
    end

end

pe_out = reshape(pe,size(t_mev_in));
ee_out = reshape(ee,size(t_mev_in));
se_out = reshape(se,size(t_mev_in));
ue_out = reshape(ue,size(t_mev_in));
yeplus_out = reshape(yeplus,size(t_mev_in));
rel_out = reshape(rel,size(t_mev_in));
rel1_out = reshape(rel1,size(t_mev_in));
rel2_out = reshape(rel2,size(t_mev_in));
approx_out = reshape(approx_out,size(t_mev_in));
    
end