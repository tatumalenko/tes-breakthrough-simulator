function ysol = Scenario1_MOL

R  = Const('R'); % J/mol.K or Pa.m^3/mol.K
R_Latm = 0.082;
PI = Const('pi');
kB = Const('k');
g  = Const('g');

% BED CHARACTERISTICS
H     = 0.0695;        % m; height of column
% eb = 0.4 + 0.05*(Sol.Dp/Sol.Dint) + 0.412*(Sol.Dp/Sol.Dint)^2; % Sol.Dp/Sol.Dint < 0.5 (Dixon 1988)
eb    = 0.39;       % --; bed void fraction = 0.375 + 0.34*Dp/Dint (Jeshar eqn)
Dint  = 0.03391;    % m; column internal diameter
Rint  = Dint/2;          % m; radius of column
Dext  = 0.0381;     % m; column external diameter
Rext  = Dext/2;          % m; outer radius of column with a thickness of 9mm
rhob  = 900;        % (1-eb)*rhop where rhop = (1-ep)*rhos

% INLET FLOW CHARACTERISTICS
T0      = 298.15;                   % K; initial temperature,% K; ambient temperature
Tginlet = T0;
Patm    = 1.1;                      % atm; operating pressure, % Pa; atmospheric pressure (outside column)
Ppa     = Patm*101325;          % Pa; operating Pressure
yw      = 0.024;                    % mol fraction; water composition
cginlet = (yw*Patm)*1000/(0.082*T0);   % mol/m^3; inlet humidity conc.
Mw      = 18.01;                    % g/mol; water molecular weight
Ma      = 28.97;                    % g/mol; air molecular weight
Y       = (yw*(Mw/Ma))/2; % kg water/kg DA; air humidity
PartialP= Ppa*Y/((Mw/Ma)+Y); % Pa; water vapor partial pressure
Vdot    = 24;                       % LPM; volumetric flow rate
vg      = Vdot*4/(PI*(Dint^2)*60000); % m/s; bulk superficial velocity
vi      = vg/eb;                    % m/s; bulk interstitial velocity

cginlet = (yw*Patm)*1000/(R_Latm*T0); % mol/m^3; inlet humidity conc.
Tamb    = T0;

% PELLET CHARACTERISTICS
Dp   = 0.0023;   % m; mean particle diameter
Rp   = Dp/2;
rp   = 1.75e-9;  % m; mean pore radius
ep   = 0.395;     % --; particle void fraction
rhos = 2000; %1970;     % kg/m3; particle density
Cps  = 836.8;    % J/kg.K; heat capacity of solid
ks   = 0.147;    % W/m.K
rhop = 1020; %(1 - ep)*rhos;
dHads= -52000;                    % J/mol; heat of adsorption

% COLUMN WALL CHARACTERISTICS
rhow = 8238;  % kg/m3; density of column wall
Cpw  = 473;   % J/kg.K; heat capacity of the wall
kw   = 17;    % W/m.K; thermal conductivity of wall
%Uo   = 2.365;

et    = eb + (1 - eb)*ep;

DM      = 0.27; % cm2/s; molecular diffusivity
Dm      = DM/(100^2); % m2/s; " " "
Nud     = 4.36; % Nusselt number in column
nuair   = 0.00001589; % m2/s; air kinetic visc at 300K;
k       = 0.4; % W/m.K; thermal conductivity of air at ~300 K
kair    = 0.0263;% W/m.K; thermal conductivity of air at 300 K???
Prair   = 0.707; % Prandlt air number at 300K

% Temperature Fitted Paramters for Temperature-Dependent Toth
a0 = 3.634*10^(-6); % mol/kg.kPa
b0 = 2.408*10^(-7); % 1/kPa
A = 6.852*10^3; % K
t0 = 0.3974; % Tooth isotherm exponent coefficient
cons = -4.199; % K

% LDF Model Parameters
Dso = 250; % m2/s
Ea = -45500; % J/mol

% ---------------------------------------------------------------------
% GRID IN AXIAL DIRECTION
% ---------------------------------------------------------------------
nz  = 5;
dz  = H/(nz - 1);
dzs = dz^2;
z   = [0:dz:H]';
%D1  = three_point_centered_uni_D2(0,H,nz);
D1  = five_point_biased_upwind_D1(z,vg);
D2  = five_point_centered_D2(z);

% ---------------------------------------------------------------------
% INDEPENDENT VARIABLE FOR ODE INTEGRATION
% ---------------------------------------------------------------------
tf    = 8000;                % s; final time
nt    = 100;
tout  = linspace(0,tf,nt)';

% ---------------------------------------------------------------------
% INITIALIZE VARIABLES
% ---------------------------------------------------------------------
cgz=zeros(nz,1); cgzz=zeros(nz,1); cgt=zeros(nz,1);
Tgz=zeros(nz,1); Tgzz=zeros(nz,1); Tgt=zeros(nz,1);
Twt=zeros(nz,1);
qt=zeros(nz,1);

% ---------------------------------------------------------------------
% SET INITIAL CONDITIONS
% ---------------------------------------------------------------------
cg0 = zeros(nz,1);
cg0(1) = cginlet;
Tg0 = zeros(nz,1)  + T0;   %Tg(:,1)=T0;
Tw0 = zeros(nz,1)  + T0;   %Tw(:,1)=T0;
q0  = zeros(nz,1);         %q(:,1) =T0;

y0  = [cg0(:); Tg0(:); Tw0(:); q0(:)];

% ---------------------------------------------------------------------
% ODE INTEGRATION OPTIONS
% ---------------------------------------------------------------------
reltol   = 1.0e-4;
abstol   = 1.0e-4;
y0_index = [1:length(y0)];
%options  = odeset('RelTol',reltol,'AbsTol',abstol,'NonNegative',y0_index);
options  = odeset('RelTol',reltol,'AbsTol',abstol);

[t, y] = ode15s(@model,tout,y0,options); %returns mat with nvar*nz*nr cols and t rows

% % call to ODE solver
% M = mass_react;
% options = odeset( 'Mass' ,M, 'MassSingular' , 'yes' , ...
%                   'RelTol' , 1e?3, 'AbsTol' , 1e?3);
% [tout,yout] = ode15s(@model,t,x,options) ;


% ---------------------------------------------------------------------
% RETURNS MATRICES THAT ARE ntxnz*nr IN SIZE
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

Conc24LPM= xlsread('Dansdata','24LPM Conc');
Temp24LPM= xlsread('Dansdata','24LPM Temp');
e.Conc.t    = Conc24LPM(:,1);
e.Conc.cgend= Conc24LPM(:,2);
e.Temp.t    = Temp24LPM(:,1);
e.Temp.Tgend= Temp24LPM(:,2) - 273.15;

figure;
subplot(1,2,1); plot(e.Conc.t,e.Conc.cgend,'o',ysol.t,ysol.cgend); ylim([0 1.2]);
subplot(1,2,2); plot(e.Temp.t,e.Temp.Tgend,'o',ysol.t,ysol.Tgend); ylim([0 80]);

    function yt = model(t,y)
        % -----------------------------------------------------------------
        % UNPACK Y (A SINGLE ROW VECTOR) INTO EACH DEPENDENT VARIABLE
        % -----------------------------------------------------------------
        cg  = y(1:nz);          %reshape(y(1        : nz   ),nz,1);
        Tg  = y(nz+1:2*nz);     %reshape(y(nz+1     : 2*nz ),nz,1);
        Tw  = y(2*nz+1:3*nz);   %reshape(y(2*nz+1   : 3*nz ),nz,1);
        q   = y(3*nz+1:4*nz);   %reshape(y(3*nz+1   : 4*nz ),nz,1);
        
        cg(cg<0) = 0;
        q(q<0)   = 0;
        %         q(1)     = 0;
        %         cg(1) = cginlet;
        %         Tg(1) = T0;
        %         q(1) = qinlet;

        Tav  = (Tg+T0)/2;
        rhog = (Ppa-(0.378*PartialP))./(287.1.*Tav); % kg/m3; humid air dens.
        
        Cpa = 1000*(28.088+(0.197*10^(-2).*Tav)+(0.48*10^(-5).*(Tav.^2))...
            -(1.965*10^(-9)*(Tav.^3)))./Ma; % J/kg.K; DA heat capacity
        Cpv = 1000.*(32.218+(0.192.*10^(-2).*Tav)+(1.055.*10^(-5).*(Tav.^2))...
            -(3.593.*10^(-9).*(Tav.^3)))./Mw; % J/kg.K; water heat capacity
        Cpg = Cpa+(Cpv.*Y); % J/kg.K; humid air heat capacity;
        
        mug = (1.827*10^(-5))*((291.15+120)/(Tg+120))*...
            ((Tg./291.15).^(3/2)); % kg/m.s; air dvisc
        
        Rep = rhog.*vg.*Dp./mug; % Reynolds particle number
        Sc  = mug./(rhog.*Dm); % Schimdt gas number
        Dz  = (Dm./eb).*(20.0+0.5.*Rep.*Sc); % m^2/s; axial dispersion coeff
        %kg  = (Dm/Dp)*(2.0+1.1*(Rep^0.6)*(Sc^0.33)); % m/s; mass t.c
        hw  = k.*Nud./Dint; % column overall heat t.c
        Rad = (1/298.15)*(298.15-T0).*g ...
            .*(Dext^3).*Prair./(nuair^2); % Rayleight air number
        Nu  = (0.60+0.387*(Rad^(1/6))*(1 ...
            +(0.559./Prair)^(9/16))^(-8/27))^2; % Nusselt number outside
        ho  = kair.*Nu./(Dext); % Overall heat t.c outside
        eta = rhog.*Cpg.*(eb+(1-eb)*ep)+...
            (1-eb)*(1-ep).*rhos.*Cps; % weighted average heat capacity
        
        b   = b0.*exp(A./Tg); % 1/kPa
        a   = a0.*exp(A./Tg); % mol/kg.kPa
        n   = t0 + cons./Tg;
        cge = q.*101.3./(((a.*R./1000.*Tg).^n) ...
            -((q.*b.*R./1000.*Tg).^n)).^(1./n); % mol/m^3; rearranged Toth
        
        kads= (15.*Dso.*exp(Ea./(8.314.*(Tg))))./(rp.^2);

        for iz=1:nz
            if(iz==1)
                %cg_bcw = -vi/Dz*(cginlet - cg(i));
                %Tg_bcw = -rhog*Cpg*vg/k*(Tginlet - Tg(i)); %changed vg to vi
                %cg_w   = cg(i+1) - 2*(dz)*cg_bcw;% cg(i) - (dz)*cg_bcw; %
                cg_w   = cginlet;
                cg_e   = cg(iz+1);             %
                Tg_w   = Tginlet; %Tg(i+1) - 2*(dz)*Tg_bcw;% Tg(i) - (dz)*Tg_bcw; %
                Tg_e   = Tg(iz+1);
                
            elseif(iz==nz)
                cg_bce = 0;
                Tg_bce = 0;
                cg_w   = cg(iz-1);
                cg_e   = cg(iz) + (dz)*cg_bce;%cg(i-1) - 2*(dz)*cg_bce; %cg(i) + (dz)*cg_bce;
                Tg_w   = Tg(iz-1);
                Tg_e   = Tg(iz) + (dz)*Tg_bce;%Tg(i-1) - 2*(dz)*Tg_bce; %Tg(i) + (dz)*Tg_bce;
            else
                cg_w = cg(iz-1);
                cg_e = cg(iz+1);
                Tg_w = Tg(iz-1);
                Tg_e = Tg(iz+1);
            end
            cgz(iz)  = (cg(iz) - cg_w)/dz;
            cgzz(iz) = (cg_e - 2.0*cg(iz) + cg_w)/dzs;
            Tgz(iz)  = (Tg(iz) - Tg_w)/dz;
            Tgzz(iz) = (Tg_e - 2.0*Tg(iz) + Tg_w)/dzs;
            
            %             % first and second order derivative
            %             cgz = D1*cg;
            %             Tgz = D1*Tg;
            %             cgzz= D2*cg;
            %             Tgzz= D2*Tg;
            
            % temporal derivatives
            %             cgt = -vg*cgz + Dz*cgzz - kads*(cg(i)-cge);
            %             Tgt = -vg*Tgz + k*Tgzz - 2/Rint*hw/eta*(Tg-Tw) - dHads*kads*eb/eta.*(cg-cge);
            %             Twt = 2*Rint*hw/((Rext^2 - Rint^2)*rhow*Cpw).*(Tg-Tw) ...
            %                 -2*Rext*ho/((Rext^2 - Rint^2)*rhow*Cpw).*(Tw-Tamb);
            %             qt  = (kads*eb/rhop)*(cg-cge);
            
            % boundary conditions at z = z01 cB1t(1) =cBin?cB1(1);
%             cgt(1) = cginlet - cg(1);
%             Tgt(1) = Tginlet - Tg(1);
%             
%             % boundary conditions at z = zL1 = z02 cB1t(n1) = cB2z(1) ? cB1z(n1);
%             cgt(nz) = -cgz(nz);
%             Tgt(nz) = -Tgz(nz);
            
            % -------------------------------------------------------------
            % PDES
            % -------------------------------------------------------------
            %dqdt = (kads*eb/rhob)*(cg(i)-cge);
            cgt(iz)=...
                Dz(iz)                      *cgzz(iz)...         Diffusive
                -vi                         *cgz(iz)...          Convective
                -kads(iz)                   *(cg(iz)-cge(iz))...     Source
                ;
            
            Tgt(iz)=...
                k/eta(iz)                      *Tgzz(iz)...         Diffusive
                -vg.*rhog(iz).*Cpg(iz)./eta(iz)            *Tgz(iz)...          Convective
                -2/Rint*hw/eta(iz)              *(Tg(iz)-Tw(iz))...   Source
                -dHads.*kads(iz).*eb./eta(iz)         *(cg(iz)-cge(iz))...     Source
                ;
            
            Twt(iz)=...
                2*Rint*hw/((Rext^2 - Rint^2)*rhow*Cpw)  *(Tg(iz)-Tw(iz))... Source
                -2*Rext*ho/((Rext^2 - Rint^2)*rhow*Cpw) *(Tw(iz)-Tamb)...  Source
                ;
            
            qt(iz)=...
                (kads(iz)*eb/rhop)*(cg(iz)-cge(iz));
            
            %             Md = spdiags([1/dzs.*ones(nz,1) -2/dzs.*ones(nz,1) 1/dzs.*ones(nz,1)],[-1 0 1],nz,nz);
            %             Mc = spdiags([1/dz.*ones(nz,1) -1/dz.*ones(nz,1)],[0 -1],nz,nz);
            %
            %             transient = [1; 1; 1; 1];
            %             diffusive = [Dz; k/eta; 0; 0];
            %             convective = [-vi; -vg*rhog*Cpg/eta; 0; 0];
            %             source = [-kads*(cg(i)-cge); ...
            %                       -2/Rint*hw/eta*(Tg(i)-Tw(i))-dHads*kads*eb/eta*(cg(i)-cge)
            %                       2*Rint*hw/((Rext^2 - Rint^2)*rhow*Cpw)*(Tg(i)-Tw(i)) - 2*Rext*ho/((Rext^2 - Rint^2)*rhow*Cpw)*(Tw(i)-Tamb); ...
            %                       (kads*eb/rhop)*(cg(i)-cge)];
        end
        % -----------------------------------------------------------------
        % 2D to 1D matrices
        % -----------------------------------------------------------------
        yt    = [cgt; Tgt; Twt; qt];
        
    end
% 
%     function M = mass_diffusion(n)
%         % Mass matrix
%         M = eye(n);
%         M(1,1) = 0;
%         M(n,n) = 0;
%         M = sparse(M);
%     end

end








