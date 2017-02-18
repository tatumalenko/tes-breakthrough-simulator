function dqdc = LI2_isotherm(c,T,k1,k2,k3,k4)
%     if c < 0
%         c = 0;
%     end
    T  = T + 273.15;
    R  = R_GAS('kPa','m3');    % atm.m3/molG.K
    
    P  = c.*R.*T; % molG/m3G*atm.m3G/molG/K*K = atm
    
%     Mv = 18.01528e-3;         % kg/mol
    
%     Hsv     = CoolProp.PropsSI('H','T',T,'Q',1,'Water');
%     Hsl     = CoolProp.PropsSI('H','T',T,'Q',0,'Water');
%     dHvap   = (Hsl - Hsv).*Mv; % J/mol
%     
%     k1 = 174.2;             % molA/kgS
%     k2 = 0.510;             % molA/kgS.K
%     k3 = 1.375;             % 1/atm
%     k4 = 408.2;             % K 
%     
    qm = (k1 - k2*T);       % molA/kgS
    b  = k3.*exp(k4./T);    % 1/atm
    q  = (qm.*b.*P)./(1 + b.*P); % molA/kgS
    
%     Qst      = R_GAS('J').*k4;         % J/mol
%     dHads    = Qst - dHvap;            % J/mol

%     theta    = q./qm.*100;
%     W        = q.*Mv.*100;              % kgA/kgS*100
    dqdp     = (qm.*b)./(1 + b.*P).^2; % molA/kgS.atm
    dqdc     = dqdp.*R.*T;             % (molA/molG)*m3G/kgS = molA/kgS*m3G/molG
    
%     p1 = -0.1348;
%     p2 = 7.578;
%     p3 = -183.9;
%     p4 = 4854;
%     dh = (p1.*W.^3 + p2.*W.^2 + p3.*W + p4).*1000.*Mv; %J/g*g/kg*kg/mol=J/mol
%     
%     LI.theta = theta;
%     LI.W     = W;
%     LI.dh    = dh;
% %     LI.dHads = dHads;
    LI.qmol  = q; 
%     LI.qkg   = q.*Mv;
    LI.dqdp  = dqdp;
    LI.dqdc  = dqdc;
end