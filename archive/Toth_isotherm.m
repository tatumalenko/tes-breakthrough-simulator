function Toth = Toth_isotherm(cp,Tp)
    % Temperature Fitted Paramters for Temperature-Dependent Toth
    R = 0.008313;           % m3.kPa/mol.K
    a0 = 3.634*10^(-6);     % mol/kg.kPa
    b0 = 2.408*10^(-7);     % 1/kPa
    A = 6.852*10^3;         % K
    t0 = 0.3974;            % Toth isotherm exponent coefficient
    cons = -4.199;          % K
    b = b0*exp(A/Tp);       % 1/kPa
    a = a0*exp(A/Tp);       % mol/kg.kPa
    n = t0 + cons/Tp;

    P    = cp*R*Tp;
    q    = a*P/(1 + (b*P)^n)^(1/n);
    dqdp = a*((b*P)^n + 1)^(-(n+1)/n);
    dqdc = R*Tp*dqdp;
    
%     if dqdc > 100
%         Toth.dqdc = 100;
%     else
    Toth.q    = q;
    Toth.dqdp = dqdp;
    Toth.dqdc = dqdc;
%     end
end