function Sol = transferCoeffs(Tg,Tw,Pars)
    Rg    = Pars.R;         % J/mol.K; Gas constant
    PI    = Pars.PI;        % --; Pi constant
    kB    = Pars.kB;        % m^2.kg/s^2.K; Boltzmann constant
    g     = Pars.g;         % m/s^2; Gravitational constant
    H     = Pars.H;         % m; bed height
    Dint  = Pars.Dint;      % m; inside column diameter
    Dext  = Pars.Dext;      % m; outside column diameter
    kw    = Pars.kw;        % W/m.K; wall thermal conductivity
    ks    = Pars.ks;
    Dp    = Pars.Dp;        % m; mean particle diameter
    rp    = Pars.rp;        % m; mean pore radius
    Tamb  = Pars.T0;        % K; ambient temperature
    Patm  = Pars.Patm;      % Pa; atmospheric pressure (outside column)
    Ptot  = Pars.Ppa;
%     RHi   = Pars.RHi;
%     RHo   = Pars.RHo;
    yw    = Pars.yw;
    vg    = Pars.vg;        % m/s; bulk interstitial velocity
    Ma    = 28.9635e-3;     % kg/mol
    Mv    = 18.01528e-3;    % kg/mol
    dipmo_v = Pars.dipmo_v; % debyes; water dipole moment
    dipmo_a = Pars.dipmo_a; % debyes; air dipole moment
    Vlmnb_v = Pars.Vlmnb_v; % cm^3/mol; water liquid molar volume @ Tnb
    Vlmnb_a = Pars.Vlmnb_a; % cm^3/mol; air liquid molar volume @ Tn
    Tnb_v = Pars.Tnb_v;     % K; water normal boiling temperature
    Tnb_a = Pars.Tnb_a;     % K; air normal boiling temperature
    eb    = Pars.eb;        % --; bed void fraction = 0.375 + 0.34*Dp/Dint (Jeshar eqn)
    ep    = Pars.ep;        % --; particle void fraction
    
    % FLOW PROPERTIES
    Tbulk = Tg;             % K; current bulk temperature
    Twall = Tw;             % K; current wall temperature
    
    % FLUID PROPERTIES @ Tav (bulk-wall avg)
    Tav     = (Twall+Tbulk)/2;
    BulkAir = propAir(Tav,Ptot,'yw',yw);
    rhog    = BulkAir.rho;
    dviscg  = BulkAir.dvisc;
    kviscg  = BulkAir.kvisc;
    kg      = BulkAir.k;
    Cpg     = BulkAir.Cp;
    %     alphag  = BulkAir.alpha;
    %     betag   = BulkAir.beta;
    Prg     = BulkAir.Pr;
    Rep     = rhog*eb*vg*Dp/dviscg;     % Reynolds particle number
        
    % BULK AND PELLET MT PROPERTIES
    tort = ep + 1.5*(1-ep);             % --; tortuosity
    Mav  = 2*(1/Ma + 1/Mv)^-1;          % kg/mol; some sort of MW average
    %----------------------------------------------------------------------
    %% CHAPMANN-ENSKOG
    %----------------------------------------------------------------------
    delta_a = 1.94e3*dipmo_a^2/(Vlmnb_a*Tnb_a); % debyes^2.mol/cm^3.K
    delta_v = 1.94e3*dipmo_v^2/(Vlmnb_v*Tnb_v);
    delta_av = (delta_a*delta_v)^0.5;
    
    % eps_akB = 78.6*kB;                        % air characteristic L-J energy, 78.6 = eps_a*kB (K)
    eps_akB = 1.18*(1 + 1.3*delta_a^2)*Tnb_a;
    % eps_vkB = 809.1*kB;                       % vapour characteristic L-J energy/kB, 809.1 = eps_a*kB (K)
    eps_vkB = 1.18*(1 + 1.3*delta_v^2)*Tnb_v;
    eps_avkB = (eps_akB*eps_vkB)^0.5;
    
    % sigma_a = 3.711e-10;                      % A; collision diameter of air
    sigma_a = (1.585*Vlmnb_a/(1+1.3*delta_a^2))^(1/3);
    % sigma_v = 2.641e-10;                      % A; collision diameter of vapour
    sigma_v = (1.585*Vlmnb_v/(1+1.3*delta_v^2))^(1/3);
    sigma_av = 10^-10*(sigma_a*sigma_v)^0.5;    % m; geometric average A->m.10^-10
    
    Tstar = Tg/eps_avkB;

    Omega = 1.06036/Tstar^0.15610+0.19300/exp(0.47635*Tstar)+1.03587/exp(1.52996*Tstar)+1.76474/exp(3.89411*Tstar); % collision integral relation by Neufield 
    OmegaP = Omega + 0.19*delta_av^2/Tstar;
    
    zg   = 1;                                   % --; bulk gas compressibility factor
    nd   = Ptot/(zg*Rg*Tg);                     % mol/m^3; number (molar) density
    Dm2   = (3/16)*(4*PI*kB*Tg/Mav)^0.5/(nd*PI*(sigma_av)^2*OmegaP); % m^2/s; molecular diffusivity
    
    %----------------------------------------------------------------------
    %% FULLER ET AL. METHOD 
    %----------------------------------------------------------------------
    vD_a = 19.7;                        % Diffusion volumes of air
    vD_v = 13.1;                        % Diffusion volumes of water
    Dm3   = 0.00143*Tg^1.75/((Ptot*10^-5)*(Mav*1000)^0.5*(vD_a^(1/3)+vD_v^(1/3))^2)/100^2; % cm^2->m^2/s; molecular diffusivity
    
    %----------------------------------------------------------------------
    %% EXPERIMENTAL VALUES FOR H2O-AIR @ 313 K. Ref: Carmichael, et al. (1955)
    %----------------------------------------------------------------------
    DmP  = 0.292*10;                    % DmP[(cm2/s).bar]*10^5/(100^2)=DmP[(m2/s).Pa]
    Dm  = DmP*(Tg/313)^1.5/Ptot;       % Dm corrected for pressure and temperature
    
    %% ----------------------------------------------------------------------
    % CALCULATIONS
    % ----------------------------------------------------------------------
    %% MASS TRANSFER CORRELATIONS
    Dk   = 2/3*rp*(8*Rg*Tbulk/(PI*Mv))^0.5; % m/s; knudsen diffusivity
    De   = (1/Dk + 1/Dm)^(-1)/tort;
    Sc   = dviscg/(rhog*Dm);            % Schimdt gas number
    Dz   = (Dm/eb)*(20 + 0.5*Rep*Sc);   % m^2/s; axial dispersion coeff Wakao (19
    % Dz = 
    Sh   = 2 + 1.1*Sc^(1/3)*Rep^(0.6);  % 3 < Re < 10^4
    kf   = (Dm/Dp)*Sh;                  % m/s; mass t.c
    Pe   = Dz/(vg*Dp);                  % 20/(Sc*Re) + 0.5;
    % Bim = kf*Rp/(3*ep*De);               % = Sh/6*Dm/(ep*De);
    
    %% BULK AND PELLET HT CORRELATIONS
    % Bih = hf*(Dp/2)/(3*ks); % = Nuf/6*kg/ks;
    Lz = kg*(7 + 0.5*Prg*Rep);
    n      = 0.28 - 0.757*log10(eb) - 0.057*log10(ks/kg);
    ke0    = kg*(ks/kg)^n;
    kz     = ke0 + 0.5*Rep*Prg*kg;        % Wein et al. (2000)
    if Rep < 4000
        Nug = 2 + 1.1*Rep^0.6*Prg^(1/3);                  % Re < 4000
    else
        Nug = Rep*Pr^(1/3)*(20.4*Rep^-0.815/(0.95*eb));   % 5000 < Re < 10300
    end
    
    hf     = (kg/Dp)*Nug;                                 % ; film htc

    %% FLUID PROPERTIES CORRELATIONS @ Tf (ambient-wall avg)
    Tf     = (Twall + Tamb)/2;
    ExtAir = propAir(Tf,Ptot,'yw',0);
    rhoo   = ExtAir.rho;
    dvisco = ExtAir.dvisc;
    kvisco = ExtAir.kvisc;
    ko     = ExtAir.k;
    Cpo    = ExtAir.Cp;
    alphao = ExtAir.alpha;
    betao  = ExtAir.beta;
    Pro    = ExtAir.Pr;  
 
    %% WALL HT PROPERTIES CORRELATIONS
    % Nui = 0.813*Rep^0.19*exp(-12*Rp/Dint); % Leva's correlation - Ruthven (1984)
    % Nui = 0.17*Rep^0.79; % --; Li & Finlayson (1977)
    Nui    = 12.5 + 0.048*Rep;          % --; De Wasch & Froment (1972)
    hi     = (kg/Dp)*Nui;               % ; col wall htc
    
    Ra     = g*betao*((Twall - Tf)/kvisco*alphao)*H^3;
    Nuo    = 0.68 + 0.67*Ra^(1/4)/(1 + (0.492/Pro)^(9/12))^(4/9);
    ho     = (ko/H)*Nuo;
    
    Uo     = (1/hi + Dint/kw*log(Dext/Dint) + Dint/(Dext*ho))^(-1);
    
    % OUTPUT RESULTS
    Sol.Psv  = BulkAir.Psv;
    Sol.rhog = rhog;
    Sol.Cpg  = Cpg;
    Sol.De   = De; %2.4e-11;
    Sol.Dz   = Dz; %5.9e-20;
    Sol.kf   = kf; %8.2e-10;
    Sol.kg   = kg;
    Sol.Lz   = Lz;
    Sol.kz   = kz;
    Sol.hf   = hf; %700;
    Sol.hi   = hi;
    Sol.ho   = ho;
    Sol.Uo   = Uo;
    
end

function sol = propAir(T,totP,str,relHum)
Rg = 8.314; % J/mol.K
tk = 273.15; % K
Ma = 28.9635e-3; % kg/mol
Mv = 18.01528e-3; % kg/mol
for i = 1:length(relHum)
    Psv = (0.7073034146 + -2.703615165e-2.*(T - tk) + 4.36088211e-3.*(T - tk).^2 + -4.662575642e-5.*(T - tk).^3 ...
               + 1.034693708e-6.*(T - tk).^4).*1000; % Pa; vap sat pressure
    if strcmpi(str,'RH')
        RH = relHum(i);
    elseif strcmpi(str,'yw')
        RH = relHum(i).*totP./Psv;
    end

    for j = 1:length(totP)
        Po = totP(j); % Pa
                
        dvisca = (-9.8601e-1 + 9.080125e-2.*T + -1.17635575e-4.*T.^2 + 1.2349703e-7.*T.^3 + -5.7971299e-11.*T.^4)*10^-6; % N.s/m2
        dviscv = sqrt(T./647.27)./(0.0181583+0.0177624.*(647.27./T)+0.0105287*(647.27./T).^2-0.0036744*(647.27./T).^3)*10^-6; % N.s/m2
        
        ka    = (-2.276501e-3 + 1.2598485e-4.*T + -1.4815235e-7.*T.^2 + 1.73550646e-10.*T.^3 + -1.066657e-13.*T.^4 + 2.47663035e-17.*T.^5); % W/m.K
        kv    = (1.761758242e1 + 5.558941059e-2.*(T-tk) + 1.663336663e-4.*(T-tk).^2).*10.^-3; % W/m.K
        
        Cpa   = (0.103409e1 + -0.284887e-3.*T + 0.7816818e-6.*T.^2 + -0.4970786e-9.*T.^3 + 0.1077024e-12.*T.^4).*1000; % J/kg.K
        Cpv   = (1.86910989 + -2.578421578e-4.*(T-tk) + 1.941058941e-5.*(T-tk).^2).*1000; % J./kg.K
        
        Psv   = (0.7073034146 + -2.703615165e-2.*(T - tk) + 4.36088211e-3.*(T - tk).^2 + -4.662575642e-5.*(T - tk).^3 ...
               + 1.034693708e-6.*(T - tk).^4).*1000; % Pa; vap sat pressure
        
        chi1  = 3.53624e-4 + 2.93228e-5.*(T-tk) + 2.61474e-7.*(T-tk).^2 + 8.57538e-9.*(T-tk).^3;
        chi2  = exp(-1.07588e1 + 6.32529e-2.*(T-tk) + -2.53591e-4.*(T-tk).^2 + 6.33784e-7.*(T-tk).^3);
        fTP   = exp(chi1.*(1-Psv./Po) + chi2.*(Psv./Po - 1));

        xv    = fTP.*RH.*Psv./Po;
        
        zv    = 1 + (0.7e-8 + -0.147184e-8.*exp(1734.29./T)).*Psv +...
               (0.104e-14 + -0.335297e-17.*exp(3645.09./T)).*Psv.^2;
        
        rho   = 1./zv.*Po./(Rg.*T).*Ma.*(1 - xv.*(1 - Mv./Ma));
        
        phiav = sqrt(2)./4.*(1+Ma./Mv).^(-1./2).*(1+(dvisca./dviscv).^0.5.*(Mv./Ma).^0.25).^2;
        phiva = sqrt(2)./4.*(1+Mv./Ma).^(-1./2).*(1+(dviscv./dvisca).^0.5.*(Ma./Mv).^0.25).^2;
        
        dvisc = (1-xv).*dvisca./((1-xv)+xv.*phiav) + xv.*dviscv./(xv+(1-xv).*phiva);
 
        k     = (1-xv).*ka./((1-xv)+xv.*phiav) + xv.*kv./(xv+(1-xv).*phiva); % W/m.K; thermal conductivity
        
        Cp    = (Cpa.*(1-xv).*Ma+Cpv.*xv.*Mv)./(Ma.*(1-xv)+Mv.*xv);
        
        alpha = k./(rho.*Cp);
        
        beta  = 1/T;
        
        Pr    = dvisc.*Cp./k;
        
        kvisc = dvisc./rho;
        
        sol.Psv     = Psv;
        sol.xv      = xv;
        sol.zv      = zv;
        sol.rho     = rho;
        sol.dvisc   = dvisc;
        sol.kvisc   = kvisc;
        sol.k       = k;
        sol.Cp      = Cp;
        sol.alpha   = alpha;
        sol.beta    = beta;
        sol.Pr      = Pr;        
        
%         sol.Psv = Psv;
%         sol.xv{i,j} = xv;
%         sol.zv{i,j} = zv;
%         sol.rho{i,j} = rho;
%         sol.dvisc{i,j} = dvisc;
%         sol.kvisc{i,j} = kvisc;
%         sol.k{i,j} = k;
%         sol.Cp{i,j} = Cp;
%         sol.alpha{i,j} = alpha;
%         sol.beta{i,j} = beta;
%         sol.Pr{i,j} = Pr;
        
       
    end
end

end

