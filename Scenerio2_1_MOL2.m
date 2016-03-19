function out = Scenerio2_1_MOL2
%%  Scenerio2 FUNCTION
%   
    T0      = 298.15;
    
    pelDIA = 0.0023;
    pelRAD = pelDIA/2;
    NI     = 4;
    dr     = pelRAD/(NI-1);
    r      = 0:dr:pelRAD;
    
    bedHEIGHT = 0.0695;
    NK   = 20;
    dz   = bedHEIGHT/(NK-1);
    z    = 0:dz:bedHEIGHT;
    
    timeFinal = 3600;
    NT        = 100;
    
    Mesh = struct('T0',T0,'pelRAD',pelRAD,'NI',NI,'dr',dr,'bedHEIGHT',bedHEIGHT,...
                  'NK',NK,'dz',dz,'r',r,'z',z);
        
    t = linspace(0,timeFinal,NT);
    cg = zeros(NK,1);
    Tg = zeros(NK,1) + T0;
    Tw = zeros(NK,1) + T0;
    cp = zeros(NI,NK);
    Tp = zeros(NI,NK) + T0;
    y  = [cg; Tg; Tw; cp(:); Tp(:)];
    
    fmodel = @(t,y)MODEL(t,y,Mesh);
    options = odeset('RelTol',1e-3,'AbsTol',1e-3,'MaxOrder',2);
    [tout,yout] = ode15s(fmodel,t,y,options);
    
    cg = reshape(yout(1:NT,1:NK),NT,NK,1);
    Tg = reshape(yout(1:NT,(NK+1):(2*NK)),NT,NK,1);
    Tw = reshape(yout(1:NT,(2*NK+1):(3*NK)),NT,NK,1);
    cp = reshape(yout(1:NT,(3*NK+1):(end-NI*NK)),NT,NI,NK);
    Tp = reshape(yout(1:NT,(end-NI*NK+1):end),NT,NI,NK);
    
    out.y = yout;
    out.t = tout;
    out.r = r;
    out.z = z;
    out.cg = cg;
    out.cp = cp;
    out.Tg = Tg;
    out.Tp = Tp;
    out.Tw = Tw;
end

function yt = MODEL(t,y,Mesh)
    
    pelRAD      = Mesh.pelRAD;
    NI          = Mesh.NI;
    dr          = Mesh.dr;
    bedHEIGHT   = Mesh.bedHEIGHT;
    NK          = Mesh.NK;
    dz          = Mesh.dz;
    r           = Mesh.r;
    T0          = Mesh.T0;
    
    cg = reshape(y(1:NK),NK,1);
    Tg = reshape(y((NK+1):(2*NK)),NK,1);
    Tw = reshape(y((2*NK+1):(3*NK)),NK,1);
    cp = reshape(y((3*NK+1):(end-NI*NK)),NI,NK);
    Tp = reshape(y((end-NI*NK+1):end),NI,NK);

    % ------------------------------------------------------------
    % MODEL PARAMETERS
    % ------------------------------------------------------------    
    PI    = pi;
    GRAV  = 9.81;        % (m/s^2)
    RGAS  = 0.08206;     % (L.atm/mol.K)
    RGAS2 = 0.008313;    % (m3.kPa/mol.K)  
       
    yw      = 0.024;        % mol fraction; water composition

    TgMAX   = 310;          % K; to allow rough approx. of TgAVE
    Patm    = 1.1;          % atm; operating pressure
    Ppa     = Patm*101325;  % Pa; operating Pressure 
    
    cgFEED  = (yw*Patm)*1000/(RGAS*T0);     % mol/m^3; inlet humidity conc.
    TgFEED  = T0; 
    ambTg   = T0;
    TgAVE   = (TgMAX + T0)/2;
    waterMW = 18.01;                        % g/mol; water molecular weight
    airMW   = 28.97;                        % g/mol; air molecular weight
    Y       = (yw*(waterMW/airMW))/2;       % kg water/kg DA; air humidity
    partialP= Ppa*Y/((waterMW/airMW)+Y);    % Pa; water vapor partial pressure
    
    % ------------------------------------------------------------
    % Bed and particles properties
    bedDIA   = 0.03391;     % m; column internal diameter 
    bedDIAO  = 0.0381;      % m; column external diameter 
    bedRAD   = bedDIA/2;    % m; radius of column 
    bedRADO  = bedDIAO/2;   % m; outer radius of column with a thickness of 9mm
    walAREA  = (bedRADO^2)-(bedRAD^2);
    bedVOID  = 0.39;        % void fraction in the column         --> = dependant on geometry of particle, if perfect spheres: ec=0.4
    pelDIA   = 2*pelRAD;
    pelVOID  = 0.395;       % void fraction in the pellets        --> = 1-rhop/rhos (ie. rhop=(1-ep)*rhos)
    pelSPHE  = 0.9;
    pelRHO   = 1020;        % kg/m^3; solid or particle density   --> = rhob/(1-ec)
    solRHO   = 2000;        % kg/m^3; pellet bulk density         --> = property of solid (skeletal density)
    bedRHO   = 1020;        % kg/m^3; solid or particle density   --> = mass/volBed where volBed=pi*rc^2*H
    walRHO   = 8238;        % kg/m^3; density of column wall
    gasRHO   = (Ppa-(0.378*partialP))/(287.1*TgAVE);  % kg/m^3; humid air dens.
    
    % ------------------------------------------------------------
    % Thermal properties
    gasCOND = 0.4;        % W/m.K; thermal conductivity of air at ~300 K 
    ambCOND = 0.0263;     % W/m.K; thermal conductivity of air at 300 K???
    gpCOND  = 0.0146;     % W/m.K; thermal conductivity of gas-solid interface
    solCP   = 836.8;    % J/kg.K; heat capacity of solid 
    walCP   = 473;      % J/kg.K; heat capacity of the wall         
    airCP   = 1000*(28.088+(0.197*10^(-2)*TgAVE)+(0.48*10^(-5)*(TgAVE^2))...
              -(1.965*10^(-9)*(TgAVE^3)))/airMW;    % J/kg.K; DA heat capacity
    vaporCP = 1000*(32.218+(0.192*10^(-2)*TgAVE)+(1.055*10^(-5)*(TgAVE^2))...
              -(3.593*10^(-9)*(TgAVE^3)))/waterMW;  % J/kg.K; water heat capacity
    gasCP   = airCP+(vaporCP*Y);                    % J/kg.K; humid air heat capacity
    gpCP = gasCP*gasRHO*pelVOID-solCP*(pelVOID-1)*solRHO;
    % ------------------------------------------------------------
    % Pressure-velocity field
    gasFLOW = 24;                      % LPM; volumetric flow rate
    gasVEL0 = gasFLOW*4/(PI*(bedDIA^2)*6e4);   % m/s; superficial velocity 
    gasSVEL = gasVEL0;                       % m/s; assumed constant
    gasIVEL = gasSVEL/bedVOID;

    % ------------------------------------------------------------
    % Mass transfer in the bulk and pellet phase
    gasDVISC = (1.827*10^(-5))*((291.15+120)/(TgAVE+120))...
                 *((TgAVE/291.15)^(3/2));          % kg/m.s; air dvisc
    TORT     = 4;
    %Dm       = 0.001858*realpow(Tp(i,j),1.5)*realsqrt(1/Mw)
    Dm       = 0.27; % cm2/s; molecular diffusivity 
    molecDIF = Dm/(100^2); % m2/s; " " "
    knudDIF  = 9700*pelRAD*realsqrt(TgAVE/waterMW);
    effDIF   = realpow(1/molecDIF+1/knudDIF,-1)/TORT;
    pelRE    = gasRHO*gasSVEL*pelDIA/gasDVISC; % Reynolds particle number
    gasSC    = gasDVISC/(gasRHO*molecDIF); % Schimdt gas number
    axialDIS = (molecDIF/bedVOID)*(20.0+0.5*pelRE*gasSC); % m^2/s; axial dispersion coeff
    gasMTC   = (molecDIF/pelDIA)*(2.0+1.1*(pelRE^0.6)*(gasSC^0.33)); % m/s; mass t.c (kf)
    % ------------------------------------------------------------
    % Energy transfer in bulk and pellet phases     
    pelHADS = 52000;                 % J/mol; heat of adsorption
    gasNU   = 4.36;                     % Nusselt number in column
    gasHTC  = gasCOND*gasNU/bedDIA;                 % column overall heat t.c 
    pelHTC  = gasCOND*gasNU/pelDIA;                 % particle overall heat t.c (Hashi=47)
    %     effCOND = gasRHO*gasCP*(bedVOID+(1-bedVOID)*pelVOID)+...
    %               (1-bedVOID)*(1-pelVOID)*solRHO*solCP; % weighted average heat capacity
    % ------------------------------------------------------------
    % Heat losses through the column wall
    ambPR    = 0.707;                  % Prandlt air number at 300K
    ambKVISC = 0.00001589;             % m2/s; air kinetic visc at 300K;
    ambRA    = (1/298.15)*(298.15-T0)*GRAV*(bedDIAO^3)...
               *ambPR/(ambKVISC^2);    % Rayleight air number
    ambNU    = (0.60+0.387*(ambRA^(1/6))*(1+(0.559/ambPR)...
               ^(9/16))^(-8/27))^2;    % Nusselt number outside
    ambHTC   = ambCOND*ambNU/(bedDIAO);% Overall heat t.c outside
    % ------------------------------------------------------------

    Pars = struct('bedDIA',bedDIA,'bedDAIO',bedDIAO,'bedRAD',bedRAD,'bedRADO',bedRADO,...
            'walAREA',walAREA,'bedVOID',bedVOID,'pedDIA',pelDIA,'pelVOID',pelVOID,...
            'pelRHO',pelRHO,'solRHO',solRHO,'bedRHO',bedRHO,'walRHO',walRHO,...
            'gasRHO',gasRHO,'gasCOND',gasCOND,'ambCOND',ambCOND,'gpCOND',gpCOND,...
            'solCP',solCP,'walCP',walCP,'gasCP',gasCP,'gasSVEL',gasSVEL,...
            'gasDVISC',gasDVISC,'effDIF',effDIF,'axialDIS',axialDIS,'gasMTC',gasMTC,...
            'pelHADS',pelHADS,'gasHTC',gasHTC,'pelHTC',pelHTC,'ambHTC',ambHTC,'gpCP',gpCP,...
            'cgFEED',cgFEED,'TgFEED',TgFEED,'ambTg',ambTg);
    
%     cg(1) = cgFEED;
%     Tg(1) = TgFEED;
%     Tw(1) = TgFEED;
%     cp(:,1) = 0;
%     Tp(:,1) = TgFEED;

%   INITIALIZE
    dcgdt = zeros(NK,1);
    dTgdt = zeros(NK,1);
    dTwdt = zeros(NK,1);
    dcpdt = zeros(NI,NK);
    dTpdt = zeros(NI,NK);
    
    
    for K = 1:NK
        dcgdt(K) = bulkMT(K,cg,cp,Mesh,Pars);
        dTgdt(K) = bulkET(K,Tg,Tp,Tw,Mesh,Pars);
        dTwdt(K) = wallET(K,Tg,Tw,Pars);
        dcpdt(:,K) = pelMT(K,cg,cp,Tp,Mesh,Pars);
        dTpdt(:,K) = pelET(K,dcpdt,cp,Tg,Tp,Mesh,Pars);
    end
    
    yt = [dcgdt; dTgdt; dTwdt; dcpdt(:); dTpdt(:)];
end

function dcgdt = bulkMT(K,cg,cp,Mesh,Pars)
    %------------------------------------------------------------------
    % BULK MASS TRANSFER EQUATIONS
    %------------------------------------------------------------------
    NI = Mesh.NI;
    NK  = Mesh.NK;
    dz = Mesh.dz;
    pelRAD = Mesh.pelRAD;
    
    axialDIS = Pars.axialDIS;
    bedVOID = Pars.bedVOID;
    gasMTC = Pars.gasMTC;
    gasSVEL = Pars.gasSVEL;
    cgFEED = Pars.cgFEED;
    
    if K == 1
        Ak = 0;
        Bk = 2*axialDIS/dz^2;
        Ck = (dz*gasSVEL - 2*axialDIS)/(2*dz^2);
        Dk = 3*gasMTC*(cg(K) - cp(NI,K))/(pelRAD*bedVOID)...
                -cgFEED*(axialDIS/dz^2 + gasSVEL/(2*dz));
    elseif K == NK
        Ak = -2*axialDIS/dz^2;
        Bk = 2*axialDIS/dz^2;
        Ck = 0;
        Dk = 3*gasMTC*(cg(K) - cp(NI,K))/(pelRAD*bedVOID);
    else
        Ak = -(2*axialDIS + dz*gasSVEL)/(2*dz^2);
        Bk = 2*axialDIS/dz^2;
        Ck = (dz*gasSVEL - 2*axialDIS)/(2*dz^2);
        Dk = 3*gasMTC*(cg(K) - cp(NI,K))/(pelRAD*bedVOID);
    end
    
    dcgdt = Ak + Bk + Ck +Dk;
    dcgdt = -dcgdt;
end

function dTgdt = bulkET(K,Tg,Tp,Tw,Mesh,Pars)
    %--------------------------------------------------------------
    % BULK HEAT TRANSFER EQUATIONS
    %--------------------------------------------------------------
    NI = Mesh.NI;
    NK  = Mesh.NK;
    dz = Mesh.dz;
    pelRAD = Mesh.pelRAD;
    
    bedRAD = Pars.bedRAD;
    bedVOID = Pars.bedVOID;
    gasRHO = Pars.gasRHO;
    gasCP = Pars.gasCP;
    gasCOND = Pars.gasCOND;
    gasSVEL = Pars.gasSVEL;
    gasHTC = Pars.gasHTC;
    TgFEED = Pars.TgFEED;
    pelHTC = Pars.pelHTC;
        
    c0 = gasRHO*gasCP*bedVOID;
    c1 = gasCOND;
    c2 = bedVOID*gasSVEL*gasRHO*gasCP;
    c3 = 3*pelHTC*(1 - bedVOID)/pelRAD;
    c4 = 2*gasHTC/bedRAD;
    
    if K == 1
        Ak = 0;
        Bk = 2*c1/(c0*dz^2);
        Ck = (c2*dz - 2*c1)/(2*c0*dz^2);
        Dk = (-2*c1*TgFEED + dz*(2*dz*(c3*Tg(K) + c4*(Tg(K) - Tw(K)))...
                - c2*TgFEED) - 2*c3*dz^2*Tp(NI,K))/(2*c0*dz^2);
    elseif K == NK
        Ak = -2*c1/(c0*dz^2);
        Bk = 2*c1/(c0*dz^2);
        Ck = 0;
        Dk = (c3*(Tg(K) - Tp(NI,K)) + c4*(Tg(K) - Tw(K)))/c0;
    else
        Ak = -(2*c1 + c2*dz)/(2*c0*dz^2);
        Bk = 2*c1/(c0*dz^2);
        Ck = (c2*dz - 2*c1)/(2*c0*dz^2);
        Dk = (c3*(Tg(K) - Tp(NI,K)) + c4*(Tg(K) - Tw(K)))/c0;
    end
    
    dTgdt = Ak + Bk + Ck +Dk;
    dTgdt = -dTgdt;
end

function dcpdt = pelMT(K,cg,cp,Tp,Mesh,Pars)
    %--------------------------------------------------------------
    % PELLET MASS TRANSFER EQUATIONS
    %--------------------------------------------------------------
    dr = Mesh.dr;
    NI = Mesh.NI;
    
    effDIF = Pars.effDIF;
    pelVOID = Pars.pelVOID;
    solRHO = Pars.solRHO;
    gasMTC = Pars.gasMTC;
    
    A = zeros(NI,1);
    B = zeros(NI,1);
    C = zeros(NI,1);
    D = zeros(NI,1);
    
    for I = 1:NI
        r    = Mesh.r(I);
        dqdcp = tothIsotherm(I,K,cp,Tp);
        alpha = effDIF*pelVOID/(pelVOID+(1-pelVOID)*solRHO*dqdcp);
        
        if I == 1
            A(I) = 0;
            B(I) = 6*alpha/dr^2;
            C(I) = -6*alpha/dr^2;
            D(I) = 0;
        elseif I == NI
            A(I) = -2*alpha/dr^2;
            B(I) = 2*alpha/dr^2;
            C(I) = 0;
            D(I) = 2*alpha*gasMTC*(dr+r)*(cg(K)-cp(NI,K))/(r*dr*effDIF);
        else
            A(I) = alpha*(dr - r)/(r*dr^2);
            B(I) = 2*alpha/dr^2;
            C(I) = -alpha*(dr + r)/(r*dr^2);
            D(I) = 0;
        end
    end
    
    dcpdt = A + B + C + D;
    dcpdt = -dcpdt;

end

function dTpdt = pelET(K,dcpdt,cp,Tg,Tp,Mesh,Pars)
    %------------------------------------------------------------------
    % PELLET HEAT TRANSFER EQUATIONS
    %------------------------------------------------------------------
    dr = Mesh.dr;
    NI = Mesh.NI;
    
    gpCOND = Pars.gpCOND;
    gpCP = Pars.gpCP;
    bedRHO = Pars.bedRHO;
    pelHADS = Pars.pelHADS;
    pelHTC = Pars.pelHTC;
    
    A = zeros(NI,1);
    B = zeros(NI,1);
    C = zeros(NI,1);
    D = zeros(NI,1);
    
    for I = 1:NI
        r = Mesh.r(I);
        dqdcp = tothIsotherm(I,K,cp,Tp);
        if I == 1
            A(I) = 0;
            B(I) = 2*gpCOND/(dr^2*gpCP);
            C(I) = -2*gpCOND/(dr^2*gpCP);
            D(I) = -bedRHO*dqdcp*dcpdt(I,K)*pelHADS/gpCP;
        elseif I == NI
            A(I) = -2*gpCOND/(dr^2*gpCP);
            B(I) = 2*gpCOND/(dr^2*gpCP);
            C(I) = 0;
            D(I) = (-bedRHO*dqdcp*dcpdt(I,K)*r*dr*pelHADS + 2*pelHTC*(dr+r)*(Tg(K)-Tp(NI,K)))/(r*dr*gpCP);
        else
            A(I) = gpCOND*(dr-r)/(r*dr^2*gpCP);
            B(I) = 2*gpCOND/(dr^2*gpCP);
            C(I) = -gpCOND*(dr+r)/(r*dr^2*gpCP);
            D(I) = -dqdcp*dcpdt(I,K)*pelHADS*bedRHO/gpCP;
        end
    end
    
    dTpdt = A + B + C +D;
    dTpdt = -dTpdt;
end

function dTwdt = wallET(K,Tg,Tw,Pars)
    %------------------------------------------------------------------
    % WALL HEAT TRANSFER EQUATIONS
    %------------------------------------------------------------------   
    
    ambHTC = Pars.ambHTC;
    gasHTC = Pars.gasHTC;
    walAREA = Pars.walAREA;
    walCP = Pars.walCP;
    walRHO = Pars.walRHO;
    ambTg  = Pars.ambTg;
    bedRADO = Pars.bedRADO;
    bedRAD = Pars.bedRAD;
    
    dTwdt = (-2*ambHTC*ambTg*bedRADO + 2*Tw(K)*(ambHTC*bedRADO + bedRAD*gasHTC)...
            -2*bedRAD*gasHTC*Tg(K))/(walAREA*walCP*walRHO);
    dTwdt = -dTwdt;
end

function dqdcp = tothIsotherm(I,K,cp,Tp)
    %------------------------------------------------------------------
    % TOTH ISOTHERM EQUATION
    %------------------------------------------------------------------
    tothA0 = 3.634*10^(-6);   % mol/kg.kPa
    tothB0 = 2.408*10^(-7);   % 1/kPa
    tothE  = 6.852*10^3;      % K
    tothN0 = 0.3974;          % exponent coefficient
    tothD  = -4.199;          % K
    
    RGAS2  = 0.008313;        % (m3.kPa/mol.K) 
    
    tothB    = tothB0*exp(tothE/Tp(I,K));                      % 1/kPa
    tothA    = tothA0*exp(tothE/Tp(I,K));                      % mol/kg.kPa
    tothN    = tothN0 + tothD/Tp(I,K);
    pelP     = RGAS2*cp(I,K)*Tp(I,K);                          % m3.kPa/(mol.K)*mol/m3*K=kPa
    dqdP     = tothA*(1+(tothB*pelP)^tothN)^(-(1+tothN)/tothN);% mol/kg.kPa
    dqdcp    = dqdP*RGAS2*Tp(I,K);
    
%     qs    = 42.70864;
%     b0    = 33.30415;
%     qr    = 0;
%     m0    = 0.6350709;
%     a     = 0;
%     c     = 297.15;
%     
%     m = m0 + a*(1 - c/Tp(I,K));
%     b = b0*exp(qr*(c/Tp(I,K) - 1));
%     dqdp  = (b*qs)/((b*cp(I,K)*RGAS2*Tp(I,K))^m + 1)^((m + 1)/m);
%     dqdcp  = dqdp*RGAS2*Tp(I,K);
end