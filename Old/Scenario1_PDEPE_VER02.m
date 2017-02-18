function ysol = Scenario1_PDEPE_VER02(varargin)
%--------------------------------------------------------------------------
if ~isempty(varargin)
    xOpt = varargin{1};
    xN   = varargin{2};
else
    xOpt = input('Enter xOpt:');
    xN = input('Enter xN:');
end
%--------------------------------------------------------------------------
global T0 cginlet PI vi Dz % declare global variables
PI = pi();
%--------------------------------------------------------------------------
Conc24LPM= xlsread('Dansdata','24LPM Conc');
Temp24LPM= xlsread('Dansdata','24LPM Temp');
e.Conc.t    = Conc24LPM(:,1);
e.Conc.cgend= Conc24LPM(:,2);
e.Conc.cgendnorm = e.Conc.cgend./max(e.Conc.cgend);
e.Temp.t    = Temp24LPM(:,1);
e.Temp.Tgend= Temp24LPM(:,2) - 273.15;
%--------------------------------------------------------------------------
tM = 300;
m = 0;      % defines geometry of space mesh (0=rectangular)
H = 0.0695; % m; height of column
t_f = 8000; % s; final time
%--------------------------------------------------------------------------
switch xOpt
    case 1
        x = linspace(0,H,xN); 
        leg = ['EquiD (x_N=',num2str(xN),')'];
    case 2
        x = [0 logspace(-1,1,xN)]./10.0.*H; 
        leg = ['Tight inlet (x_N=',num2str(xN),')'];
    case 3
        x = f_tightExit(xN,H);
        leg = ['Tight exit(x_N=',num2str(xN),')'];
    case 4
        dx = min(diff(f_tightExit(xN,H)));
        x = linspace(0,H,H/dx + 1); 
        xN = length(x);
        leg = ['Tight equiD (x_N=',num2str(xN),')'];
    case 5
        dx = max(diff(f_tightExit(xN,H)));
        x = linspace(0,H,H/dx + 1);
        xN = length(x);
        leg = ['Loose equiD (x_N=',num2str(xN),')'];
end
%--------------------------------------------------------------------------
t = linspace(0,t_f,tM);
options = odeset('AbsTol',1e-6,'RelTol',1e-6);
%--------------------------------------------------------------------------
R1 = 0.008313; % m3.kPa/mol.K
R = 0.082; % L.atm/mol.K; ideal gas law constant
T0 = 298.15; % K; initial/surrounding/inlet temperature
Tamb = T0;
Patm = 1.1; % atm; operating pressure
Ppa = Patm*101325; % Pa; operating Pressure
yw = 0.024; % mol fraction; water composition
Mw = 18.01; % g/mol; water molecular weight
Ma = 28.97; % g/mol; air molecular weight
cginlet = (yw*Patm)*1000/(R*T0); % mol/m^3; inlet humidity conc.
ec = 0.39; % void fraction in the column
ep = 0.395; % void fraction in the pellets
rhop = 1020; % kg/m^3; solid or particle density
rhos = 2000; % kg/m3; pellet bulk density
rhow = 8238; % kg/m3; density of column wall
nuair = 0.00001589; % m2/s; air kinetic visc at 300K;
Dp = 0.0023; % m; particle diameter
Dc = 0.03391; % m; column internal diameter
Dco = 0.0381; % m; column external diameter
k = 0.4; % W/m.K; thermal conductivity of air at ~300 K
kair = 0.0263;% W/m.K; thermal conductivity of air at 300 K???
Cps = 836.8; % J/kg.K; heat capacity of solid
ks  = 0.147;    % W/m.K
Cpw = 473; % J/kg.K; heat capacity of the wall
kw = 17;    % W/m.K; thermal conductivity of wall
g  = 9.81; % m/s^2; gravity
rp = Dp/2; % m; radius of pellet
rc = Dc/2; % m; radius of column
ro = Dco/2; % m; outer radius of column with a thickness of 9mm
dHads = -54000; % J/mol; heat of adsorption
Vdot = 24; % LPM; volumetric flow rate
vg0 = Vdot*4/(PI*(Dc^2)*60000); % m/s; superficial velocity
DM = 0.27; % cm2/s; molecular diffusivity
Dm = DM/(100^2); % m2/s; " " "
% Temperature Fitted Paramters for Temperature-Dependent Toth
a0 = 3.634*10^(-6); % mol/kg.kPa
b0 = 2.408*10^(-7); % 1/kPa
A = 6.852*10^3; % K
t0 = 0.3974; % Tooth isotherm exponent coefficient
cons = -4.199; % K
% LDF Model Parameters
% Dso = 250; % m2/s
% Ea = -45500; % J/mol
Dso = 150; % m2/s
Ea = -42500; % J/mol
Bihmax = 0;
Bihmin = 0;
%----------------------------------------------------------------------
try
    sol = pdepe(m,@pdefun,@pdeic,@pdebc,x,t,options);
catch ME
    sol = pdepe(m,@pdefun,@pdeic,@pdebc,x,t,options);
end
%----------------------------------------------------------------------
ysol.t  = t;
ysol.Bihmax = Bihmax;
ysol.Bihmin = Bihmin;
ysol.cg = sol(:,:,1);        % the 2D array for conc of water in gas in mol/m3
ysol.Tg = sol(:,:,2)-273.15; % the 2D array for bulk gas temp in C
ysol.Tw = sol(:,:,3)-273.15; % the 2D array for inner wall temp in C
ysol.q  = sol(:,:,4);        % the 2D array for solid's water capacity in kg/kg
ysol.cginlet    = cginlet;
ysol.cgend      = ysol.cg(:,end);
ysol.cgendnorm  = ysol.cg(:,end)./cginlet;
ysol.Tgend      = ysol.Tg(:,end);
%----------------------------------------------------------------------
% figure;
% sp1 = subplot(1,2,1); plot(e.Conc.t,e.Conc.cgendnorm,'*',ysol.t,ysol.cgendnorm,'LineWidth',3); ylim([-0.1 1]); plot(x./H.*t_f,-0.1.*ones(length(x)),'mo');
% sp2 = subplot(1,2,2); plot(e.Temp.t,e.Temp.Tgend,'*',ysol.t,ysol.Tgend,'LineWidth',3); ylim([10 80]); plot(x./H.*t_f,10.*ones(length(x)),'mo');
% title(sp1, ['^C/_C_0: ' leg], 'FontWeight','Bold');
% title(sp2, ['T_g: ' leg], 'FontWeight','Bold');

% dockTile
%----------------------------------------------------------------------
    function [c,f,s] = pdefun(x,t,u,DuDx)
    Tav = (u(3)+u(2))/2;
    Y = (yw*(Mw/Ma))/2; % kg water/kg DA; air humidity
    PartialP = Ppa*Y/((Mw/Ma)+Y); % Pa; water vapor partial pressure
    rhog = (Ppa-(0.378*PartialP))/(287.1*Tav); % kg/m3; humid air dens.

    Cpa = 1000*(28.088+(0.197*10^(-2)*Tav)+(0.48*10^(-5)*(Tav^2))-(1.965*10^(-9)*(Tav^3)))/Ma; % J/kg.K; DA heat capacity
    Cpv = 1000*(32.218+(0.192*10^(-2)*Tav)+(1.055*10^(-5)*(Tav^2))-(3.593*10^(-9)*(Tav^3)))/Mw; % J/kg.K; water heat capacity
    Cpg = Cpa+(Cpv*Y); % J/kg.K; humid air heat capacity;
    mug = (1.827*10^(-5))*((291.15+120)/(u(2)+120))*((u(2)/291.15)^(3/2)); % kg/m.s; air dvisc

    vg = vg0; % m/s; assumed constant throughout from inlet
    vi = vg/ec;
    
    Prair = 0.707; % Prandlt air number at 300K
    Pr = Prair;
    kheat = 0.03;
    
    Rep = rhog*vg*Dp/mug; % Reynolds particle number
    Rec = rhog*vg*Dc/mug;
    
    Sc = mug/(rhog*Dm); % Schimdt gas number
    Dz = (Dm/ec)*(20.0+0.5*Rep*Sc); % m^2/s; axial dispersion coeff
    %kg = (Dm/Dp)*(2.0+1.1*(Rep^0.6)*(Sc^0.33)); % m/s; mass t.c
    
    Nud = 4.36; % Nusselt number in column
    Nup = 2.0 + 1.1*Pr^0.33*Rep^0.6;
    
    hfd = k*Nud/Dc; % column overall heat t.c
    hfp = kheat*Nup/Dp; % particle overall heat t.c
    
    %Pr = Cpg*mug/k; ~0.7 kconv=0.12*Re*Pr; Rep=70.9
    
    Rad = (1/298.15)*abs(u(3)-Tamb)*g*(Dco^3)*Prair/(nuair^2); % Rayleight air number
    Nu = (0.60+0.387*(Rad^(1/6))*(1+(0.559/Prair)^(9/16))^(-8/27))^2; % Nusselt number outside
    ho = kair*Nu/(Dco); % Overall heat t.c outside
    eta = rhog*Cpg*(ec+(1-ec)*ep)+(1-ec)*(1-ep)*rhos*Cps; % weighted average heat capacity
    
    b = b0*exp(A/u(2)); % 1/kPa
    a = a0*exp(A/u(2)); % mol/kg.kPa
    n = t0 + cons/u(2);
    cpe = 101.3*u(4)/(((a*R1*u(2))^n)-((u(4)*b*R1*u(2))^n))^(1/n); % mol/m^3; rearranged Toth

    kads = (15*Dso*exp(Ea/(8.314*(u(2)))))/(rp^2);

    c = [ec; eta; ((ro^2)-(rc^2))*rhow*Cpw; 1];
    
    ks = 3.3;
    ke = ks*(1-ec)/(1+0.5*ec);
    Bih    = hfd*Dp/ke;
    
    if Bih<Bihmin || Bihmin==0
        Bihmin = Bih;
    end
    if Bih>Bihmax || Bihmax==0
        Bihmax = Bih;
    end
    
    %     f = [ec*Dz; k; 0; 0].*DuDx + [-vg; -vg*rhog*Cpg; 0; 0];
    %
    %     s = [-(kads*ec)*(u(1)-cpe); ...
    %          -(2*hfd/rc)*(u(2)-u(3)) - dHads*(kads*ec)*(u(1)-cpe); ...
    %          2*rc*hfd*(u(2)-u(3)) - 2*ro*ho*(u(3)-T0); ...
    %         (kads*ec/rhop)*(u(1)-cpe)];
    %
    f = [ec*Dz; k; 0; 0].*DuDx;

    s = [-vg; -vg*rhog*Cpg; 0; 0].*DuDx + ...
        [-(kads*ec)*(u(1)-cpe); ...
         -(2*hfd/rc)*(u(2)-u(3)) - dHads*(kads*ec)*(u(1)-cpe); ...
         2*rc*hfd*(u(2)-u(3)) - 2*ro*ho*(u(3)-T0); ...
        (kads*ec/rhop)*(u(1)-cpe)];
    end

    function u0 = pdeic(x)
        u0 = [0;T0;T0;0];
    end

    function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,t)
        % p(x,t,u)+q(x,t).*f(x,t,u,Du/Dx)=0 where the l and r represent the
        % left BCs with form: pl(0,t,u) + ql(0,t).*f(0,t,u,Du/Dx)=0
        pl = [vi/Dz*(cginlet-ul(1)); %ul(1)-cginlet; ... vi/Dz*(cginlet-ul(1));
            ul(2)-T0; ...
            ul(3)-T0; ...
            ul(4)];
        ql = [1;0;0;0];

        % right BCs with form: pr(H,t,u)+qr(H,t).*f(H,t,u,Du/Dx)=0
        pr = [0;0;0;0];
        qr = [1;1;1;1];
    end
end

function x = f_tightExit(xN,H)
    x = zeros(1,xN); 
    x(1) = 0;
    for i=2:xN
        x(i) = x(i-1) + 1/i*H; 
    end
    x = x./max(x).*H;
    x = [x(1:end-1)./max(x).*H H];
end
