function DA = DA_isotherm(c,T,Ptot)
%     if c <= 0
%         c = 0.00001;
%     end
%     assert(isreal(c)&&c>0&&isreal(T)&&T>275);
    Rg      = 8.314; % J/mol.K
    Mv      = 18.01528e-3; % kg/mol
    
    Pv      = c*Rg*T;
    Psv     = CoolProp.PropsSI('P','T',T,'Q',0,'Water');
    lnPr    = log(Psv/Pv);
    
    Hsv     = CoolProp.PropsSI('H','T',T,'Q',1,'Water');
    Hsl     = CoolProp.PropsSI('H','T',T,'Q',0,'Water');
    
    DHvap   = Hsl - Hsv;
    % DHvap1 = -7.1*T^2 + 2620.4*T + 2.2866e6; % J/kg
    
    beta    = CoolProp.PropsSI('isobaric_expansion_coefficient','P',Ptot,'T',T,'Water');
    beta20  = CoolProp.PropsSI('isobaric_expansion_coefficient','P',Ptot,'T',273.15+20,'Water');
    rho20   = CoolProp.PropsSI('Dmass','P',Ptot,'T',273.15+20,'Water');
    rho_ads = rho20/(1 + beta20*(T - 293.15));
    
    Wo      = 341.03;                % ml/g -> m3/kg;
    E       = 1192.3e3;              % J/kg;
    E2      = E*Mv;                  % J/mol;
    n       = 1.55;                  % --;
    RTE     = Rg*T/E2;
    
    A       = Rg/Mv*T*lnPr; 
    
    W       = Wo*exp(-(A/E)^n);      % m3/kg
    theta   = W/Wo;
    qs      = Wo*rho_ads/Mv;         % mol/kg
    q       = W*rho_ads/Mv;          % mol/kg
    
    dqdP_T  = n*qs*exp(-RTE^n*lnPr^n)*RTE^n*lnPr^n/(Pv*lnPr);
    dqdc    = Rg*T*dqdP_T;           % (mol_v/kg_s)/(mol_v/m3_v) = m3_v/kg_s
    
    DH      = DHvap + A;
    DH2     = DHvap + E*(log(1/theta))^(1/n) + E*beta*T/n*(log(1/theta))^-(n-1)/n; 

    DA.dqdc = dqdc;
    DA.DH   = DH;
    DA.DH2  = DH2;
end