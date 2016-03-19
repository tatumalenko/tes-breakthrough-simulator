function y = LF_isotherm(P,T)
    T  = T + 273.15;
    R  = 8.314;
    n1 = -0.3615;
    n2 = 274.23; % K
    qmax = 18; % mol/kg
    b0   = 0.000308; % 1/Pa^n
    dE   = 18016; % J/mol
    
    n = n1 + n2./T;
    b = b0.*exp(dE./(R.*T));
    q = qmax.*b.*P.^n./(1 + b.*P.^n);
    dqdP = qmax.*b.*n.*P.^(n-1)./(1 + b.*P.^n).^2;
    
    y.q = q;
    y.dqdP = dqdP;
end