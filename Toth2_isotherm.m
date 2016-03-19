function Toth = Toth2_isotherm(q,P,T)
    Mv = 0.018;
    % Temperature Fitted Paramters for Temperature-Dependent Toth
    Ru = R_GAS('m3','kPa');   % m3.kPa/mol.K
    a0 = 3.634.*10^(-6);     % mol/kg.kPa
    b0 = 2.408.*10^(-7);     % 1/kPa
    A = 6.852.*10^3;         % K
    t0 = 0.3974;            % Toth isotherm exponent coefficient
    cons = -4.199;          % K
    b = b0.*exp(A./T);       % 1/kPa
    a = a0.*exp(A./T);       % mol/kg.kPa
    n = t0 + cons./T;

    q = q./Mv;
    c_eq = q./((a.*Ru.*T).^n - (q.*b.*Ru.*T).^n).^(1./n);
%     q_eq = a.*P./(1 + (b.*P).^n).^(1/n);   % molA/kgS
%     dqdP = a.*((b.*P).^n + 1).^(-(n+1)./n);% molA/kgS.kPa
%     dqdc = R.*T.*dqdP;                     % molA/molG*m3G/kgS
    
    Toth.c_eq = c_eq;
%     Toth.q_eq = q_eq;
%     Toth.dqdP = dqdP;
%     Toth.dqdc = dqdc;

end