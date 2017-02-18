function ysol = scenario1_mol
% ---------------------------------------------------------------------
% DECLARE TEMPERATURE INDEPENDENT PARAMETERS
% ---------------------------------------------------------------------
PI = pi();          % PI constant (3.1416..)
g  = 9.81;          % m/s2; gravitational constant
R  = 8.3144621;     % J/mol.K or Pa.m^3/mol.K
R_Latm = 0.082;     % L.atm/mol.K; ideal gas law constant

% BED CHARACTERISTICS
H     = 0.0695;     % m; height of column
%eb = 0.4 + 0.05*(Sol.Dp/Sol.Dint) + 0.412*(Sol.Dp/Sol.Dint)^2; % Sol.Dp/Sol.Dint < 0.5 (Dixon 1988)
eb    = 0.39;       % --; bed void fraction = 0.375 + 0.34*Dp/Dint (Jeshar eqn)
Dint  = 0.03391;    % m; column internal diameter
Rint  = Dint/2;     % m; radius of column
Dext  = 0.0381;     % m; column external diameter
Rext  = Dext/2;     % m; outer radius of column with a thickness of 9mm
rhob  = 900;        % (1-eb)*rhop where rhop = (1-ep)*rhos

% INLET FLOW CHARACTERISTICS
T0      = 298.15;                       % K; initial temperature,% K; ambient temperature
Tginlet = T0;
Patm    = 1.1;                          % atm; operating pressure, % Pa; atmospheric pressure (outside column)
Ppa     = Patm*101325;                  % Pa; operating Pressure
yw      = 0.024;                        % mol fraction; water composition
Mw      = 18.01;                        % g/mol; water molecular weight
Ma      = 28.97;                        % g/mol; air molecular weight
Y       = (yw*(Mw/Ma))/2;               % kg water/kg DA; air humidity
PartialP= Ppa*Y/((Mw/Ma)+Y);            % Pa; water vapor partial pressure
Vdot    = 24;                           % LPM; volumetric flow rate
vg      = Vdot*4/(PI*(Dint^2)*60000);   % m/s; bulk superficial velocity
vi      = vg/eb;                        % m/s; bulk interstitial velocity

cginlet = (yw*Patm)*1000/(R_Latm*T0);   % mol/m^3; inlet humidity conc.
Tamb    = T0;

% PELLET CHARACTERISTICS
Dp   = 0.0023;  % m; mean particle diameter
Rp   = Dp/2;
rp   = Rp;      % m; mean pore radius
ep   = 0.395;   % --; particle void fraction
rhos = 2000;    % kg/m3; particle density
Cps  = 836.8;   % J/kg.K; heat capacity of solid
ks   = 0.147;   % W/m.K
rhop = 1020;    % (1 - ep)*rhos;
dHads= -54000;  % J/mol; heat of adsorption

% COLUMN WALL CHARACTERISTICS
rhow = 8238;  % kg/m3; density of column wall
Cpw  = 473;   % J/kg.K; heat capacity of the wall
kw   = 17;    % W/m.K; thermal conductivity of wall
%Uo   = 2.365;

et    = eb + (1 - eb)*ep;

DM      = 0.27;         % cm2/s; molecular diffusivity
Dm      = DM/(100^2);   % m2/s; " " "
Nud     = 4.36;         % Nusselt number in column
nuair   = 0.00001589;   % m2/s; air kinetic visc at 300K;
k       = 0.4;          % W/m.K; thermal conductivity of air at ~300 K
kair    = 0.0263;       % W/m.K; thermal conductivity of air at 300 K???
Prair   = 0.707;        % Prandlt air number at 300K

% TEMPERATURE FITTED PARAMETERS FOR TEMPERATURE-DEPENDENT TOTH
a0 = 3.634*10^(-6);     % mol/kg.kPa
b0 = 2.408*10^(-7);     % 1/kPa
A = 6.852*10^3;         % K
t0 = 0.3974;            % Toth isotherm exponent coefficient
cons = -4.199;          % K

% LDF MODEL PARAMETERS
Dso = 145;      % m2/s
Ea = -45800;    % J/mol

% ---------------------------------------------------------------------
% GRID IN AXIAL DIRECTION
% ---------------------------------------------------------------------
nz  = 50;                   % number of axial discretization nodes
dz  = H/(nz - 1);           % axial discretization step size (first deriv.)
dzs = dz^2;                 % second deriv.
z   = [0:dz:H]';            % axial discretization node array

% ---------------------------------------------------------------------
% INDEPENDENT VARIABLE FOR ODE INTEGRATION
% ---------------------------------------------------------------------
tf    = 8000;                % s; final time (end of exp. trials)
nt    = 75;                  % number of time discretization nodes
tout  = linspace(0,tf,nt)';  % time discretization node array

% ---------------------------------------------------------------------
% INITIALIZE VARIABLES
% ---------------------------------------------------------------------
cgz=zeros(nz,1); cgzz=zeros(nz,1); cgt=zeros(nz,1); % spatial and temporal cg derivative arrays
Tgz=zeros(nz,1); Tgzz=zeros(nz,1); Tgt=zeros(nz,1); % spatial and temporal Tg derivative arrays
Twt=zeros(nz,1);                                    % temporal Tw derivative array
qt=zeros(nz,1);                                     % temporal q derivative array
    
% ---------------------------------------------------------------------
% SET INITIAL CONDITIONS
% ---------------------------------------------------------------------
cg0 = zeros(nz,1);
Tg0 = zeros(nz,1)  + T0;  
Tw0 = zeros(nz,1)  + T0;   
q0  = zeros(nz,1);         

y0  = [cg0(:); Tg0(:); Tw0(:); q0(:)];

% ---------------------------------------------------------------------
% ODE INTEGRATION OPTIONS
% ---------------------------------------------------------------------
reltol   = 1.0e-6;
abstol   = 1.0e-6;
y0_index = [1:length(y0)];

options  = odeset('RelTol',reltol,'AbsTol',abstol,'NonNegative',y0_index);
[t, y] = ode15s(@model,tout,y0,options);    %returns matrix with nvar*nz cols and t rows

% ---------------------------------------------------------------------
% RETURNS MATRICES THAT ARE ntxnz IN SIZE
% ---------------------------------------------------------------------
ysol.cg  = reshape(y(1:nt,1      : nz   ),nt,nz,1);
ysol.Tg  = reshape(y(1:nt,nz+1   : 2*nz ),nt,nz,1) - 273.15;
ysol.Tgend = ysol.Tg(:,end);
ysol.Tw  = reshape(y(1:nt,2*nz+1 : 3*nz ),nt,nz,1);
ysol.q   = reshape(y(1:nt,3*nz+1 : 4*nz ),nt,nz,1);

ysol.t       = t;
ysol.y       = y;
ysol.z       = z;
ysol.cginlet = cginlet;
ysol.cgend   = ysol.cg(:,end);
ysol.cgendnorm   = ysol.cg(:,end)./cginlet;

% ---------------------------------------------------------------------
% READ AND STORE EXPERIMENTAL DATA
% ---------------------------------------------------------------------
Conc24LPM= xlsread('../test/dansdata.xlsx','24 LPM Conc');
Temp24LPM= xlsread('../test/dansdata.xlsx','24 LPM Temp');
e.Conc.t    = Conc24LPM(:,1); 
e.Conc.cgend= Conc24LPM(:,2);
e.Conc.cgendnorm = e.Conc.cgend./cginlet;
e.Temp.t    = Temp24LPM(:,1);
e.Temp.Tgend= Temp24LPM(:,2) - 273.15;

% ---------------------------------------------------------------------
% PLOT THE cg AND Tg BREAKTHROUGH CURVES AND COMPARE TO EXP. DATA
% ---------------------------------------------------------------------
figure;
subplot(1,2,1); plot(e.Conc.t,e.Conc.cgendnorm,'o',ysol.t,ysol.cgendnorm); ylim([0 1.2]);
subplot(1,2,2); plot(e.Temp.t,e.Temp.Tgend,'o',ysol.t,ysol.Tgend); ylim([0 80]);

    % ---------------------------------------------------------------------
    % DIFFERENTIAL (TEMPORAL) EQUATION ARRAY FUNCTION:
    % [dcg/dt = f(dcg/dz,dcg2/dz2,..), dTg/dt= f(dTg/dz,dTg2/dz2,..), ..]
    % ---------------------------------------------------------------------
    function yt = model(t,y)
        % Differential function with temperature changing parameters
        % recalculated at each tempo-spatial nodes 
        cg  = y(1:nz);          %reshape(y(1        : nz   ),nz,1);
        Tg  = y(nz+1:2*nz);     %reshape(y(nz+1     : 2*nz ),nz,1);
        Tw  = y(2*nz+1:3*nz);   %reshape(y(2*nz+1   : 3*nz ),nz,1);
        q   = y(3*nz+1:4*nz);   %reshape(y(3*nz+1   : 4*nz ),nz,1);
        
        cg(cg<0) = 0;
        q(q<0)   = 0;
        
        Tav  = (Tw+Tg)/2;
        rhog = (Ppa-(0.378*PartialP))./(287.1.*Tav); % kg/m3; humid air dens.
        
        Cpa = 1000*(28.088+(0.197*10^(-2).*Tav)+(0.48*10^(-5).*(Tav.^2))...
            -(1.965*10^(-9)*(Tav.^3)))./Ma;         % J/kg.K; DA heat capacity
        Cpv = 1000.*(32.218+(0.192.*10^(-2).*Tav)+(1.055.*10^(-5).*(Tav.^2))...
            -(3.593.*10^(-9).*(Tav.^3)))./Mw;       % J/kg.K; water heat capacity
        Cpg = Cpa+(Cpv.*Y);                         % J/kg.K; humid air heat capacity;
        
        mug = (1.827*10^(-5)).*((291.15+120)./(Tg+120)).*((Tg./291.15).^(3/2)); % kg/m.s; air dvisc
        
        Rep = rhog.*vg.*Dp./mug;                                    % Reynolds particle number
        Sc  = mug./(rhog.*Dm);                                      % Schimdt gas number
        Dz  = (Dm./eb).*(20.0+0.5.*Rep.*Sc);                        % m2/s; axial dispersion coeff
        % kg  = (Dm/Dp)*(2.0+1.1*(Rep^0.6)*(Sc^0.33)); % m/s; mass t.c
        hw  = k.*Nud./Dint;                                         % column overall heat t.c
        Rad = (1/T0).*abs(Tw-Tg).*g.*(Dext^3).*Prair./(nuair^2);    % Rayleight air number
        
        Nu  = (0.60+0.387.*(Rad.^(1/6))./(1+(0.559./Prair).^(9/16)).^(8/27)).^2; % Nusselt number outside
        ho  = kair.*Nu./(Dext);                                                  % Overall heat t.c outside
        eta = rhog.*Cpg.*(eb+(1-eb)*ep)+(1-eb)*(1-ep).*rhos.*Cps;                % weighted average heat capacity
        
        b   = b0.*exp(A./Tg);                       % 1/kPa
        a   = a0.*exp(A./Tg);                       % mol/kg.kPa
        n   = t0 + cons./Tg;
        cge = q.*101.3./(((a.*R./1000.*Tg).^n) ...
                -((q.*b.*R./1000.*Tg).^n)).^(1./n); % mol/m^3; rearranged Toth
        
        kads= (15.*Dso.*exp(Ea./(R.*(Tg))))./((rp).^2);
       
        % -------------------------------------------------------------
        % BOUNDARY CONDITIONS
        % -------------------------------------------------------------
        cgz_R = 0;
        cg_R  = cg(end-1) + (dz)*cgz_R;
        cg_L  = cginlet;
        
        Tgz_R = 0;
        Tg_R  = Tg(end-1) + (dz)*Tgz_R;
        Tg_L  = Tginlet;
        
        cgz(1)        = (cg(1)     - cg_L)./dz;
        cgz(2:end)    = (cg(2:end) - cg(1:end-1))./dz;
        
        cgzz(1)       = (cg_L      - 2.0*cg(1)        + cg(2))./dzs;
        cgzz(2:end-1) = (cg(3:end) - 2.0.*cg(2:end-1) + cg(1:end-2))./dzs;
        cgzz(end)     = (cg(end-1) - 2.0*cg(end)      + cg_R)./dzs;
        
        Tgz(1)        = (Tg(1)     - Tg_L)./dz;
        Tgz(2:end)    = (Tg(2:end) - Tg(1:end-1))./dz;
        
        Tgzz(1)       = (Tg_L      - 2.0*Tg(1)        + Tg(2))./dzs;
        Tgzz(2:end-1) = (Tg(3:end) - 2.0*Tg(2:end-1)  + Tg(1:end-2))./dzs;
        Tgzz(end)     = (Tg(end-1) - 2.0*Tg(end)      + Tg_R)./dzs;
        
        % -------------------------------------------------------------
        % PDES
        % -------------------------------------------------------------
        cgt=...
            Dz                          .*cgzz...         Diffusive
            -vi                         .*cgz...          Convective
            -kads                       .*(cg-cge)...     Source
            ;
        
        Tgt=...
            k./eta                       .*Tgzz...         Diffusive
            -vg.*rhog.*Cpg./eta          .*Tgz...          Convective
            -2./Rint.*hw./eta            .*(Tg-Tw)...      Source
            -dHads.*kads.*eb./eta        .*(cg-cge)...     Source
            ;
        
        Twt=...
            2.*Rint.*hw./((Rext^2 - Rint^2).*rhow.*Cpw)  .*(Tg-Tw)...   Source
            -2.*Rext.*ho./((Rext^2 - Rint^2).*rhow.*Cpw) .*(Tw-Tamb)... Source
            ;
        
        qt=...
            (kads.*eb./rhop).*(cg-cge);
        
        yt    = [cgt; Tgt; Twt; qt];
        
    end

end








