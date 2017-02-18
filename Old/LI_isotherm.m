function LI = LI_isotherm(c,T)
    if c < 0
        c = 0;
    end
    R  = R_GAS('kPa','m3');   % kPa.m3/molG.K
    P  = c.*R.*T;             % molG/m3G*kPa.m3G/molG/K*K = kPa
        
    k1 = 24.88;            % molA/kgS
    k2 = 0.0422;           % molA/kgS.K
    k3 = 1.3550e-05;       % 1/kPa
    k4 = 4865;             % K 
    
    qm = (k1 - k2*T);      % molA/kgS
    b  = k3.*exp(k4./T);   % 1/kPa
    dqdp     = (qm.*b)./(1 + b.*P).^2;   % molA/kgS.kPa
    dqdc     = dqdp.*R.*T;               % (molA/molG)*m3G/kgS = molA/kgS*m3G/molG 
     
    for sz = 1:length(dqdc)
        if dqdc(sz) > 100
            dqdc(sz) = 100;
        elseif dqdc(sz) < 1
            dqdc(sz) = 1;
        end
    end
    
    Mv = 18.01528e-3;         % kg/mol
    
%     Hsv     = CoolProp.PropsSI('H','T',T,'Q',1,'Water');
%     Hsl     = CoolProp.PropsSI('H','T',T,'Q',0,'Water');
%     dHvap   = (Hsl - Hsv).*Mv; % J/mol
    
%     Qst      = R_GAS('J').*k4;         % J/mol
%     dHads    = Qst - dHvap;            % J/mol
    q        = (qm.*b.*P)./(1 + b.*P);   % molA/kgS
    theta    = q./qm.*100;
    W        = q.*Mv.*100;               % kgA/kgS*100


    
    p1 = -0.1348;
    p2 = 7.578;
    p3 = -183.9;
    p4 = 4854;
    dh = (p1.*W.^3 + p2.*W.^2 + p3.*W + p4).*1000.*Mv; %J/g*g/kg*kg/mol=J/mol
    
    LI.theta = theta;
    LI.W     = W;
    LI.dh    = dh;
%     LI.dHads = dHads;
    LI.qmol  = q; 
    LI.qkg   = q.*Mv;
    LI.dqdp  = dqdp;
    LI.dqdc  = dqdc;
end