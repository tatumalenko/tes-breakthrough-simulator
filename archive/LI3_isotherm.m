function LI = LI3_isotherm(P,T)
    Rg = 8.314;              % J/mol.K
    Mv = 18.01528e-3;        % kg/mol
    T  = T + 273.15;
    
    %P  = c.*Rg.*T;
    k1 = 174.2;              % mol/kg
    k2 = 0.510;              % mol/kg.K
    k3 = 1.375/101325;       % 1/Pa; 1/atm -> 1/Pa (1/101325)
    %k4 = 408.2;             % K
    k4 = 310.2;              % K
    
    qm = (k1 - k2*T)./Mv;     % mol/kg -> kg_v/kg_s
    %qm = 30.*(k1 - k2.*T);   % mol/kg -> kg_v/kg_s
    b  = k3.*exp(k4./T);     % 1/Pa
    
    q        = (b.*qm.*P)./(b.*P + 1);
    dqdP     = (b.*qm)./(b.*P + 1).^2;
    LI.q     = q;
    LI.dqdc  = Rg.*T.*dqdP;
end