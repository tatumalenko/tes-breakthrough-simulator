function y = Scenerio2(Pars,varargin)
    if nargin > 1
        TinC = varargin{1};
        TinK = TinC + 273.15;
        c    = varargin{2};
        Ptot = 1.1*101325;
        for i = 1:length(TinK)
            T       = TinK(i);
            DA(i)   = DA_isotherm(c,T,Ptot);
            LI(i)   = LI_isotherm(c,T);
            TO(i)   = Toth_isotherm(c,T);
        end
        
        figure; 
        subplot(2,2,1);
        plot(TinC,dqdP);
        title('dqdP');
        subplot(2,2,2);
        plot(TinC,b);
        title('b');
        subplot(2,2,3);
        plot(TinC,dqdc);
        title('dqdc');
        subplot(2,2,4);
        plot(TinC,qm);
        title('qm');
    else
        y = Scenerio2(Pars);
    end
end

function ysol = Scenerio2(Pars)
    Rg    = Pars.R;         % J/mol.K; Gas constant
    PI    = Pars.PI;        % --; Pi constant
    kB    = Pars.kB;        % m^2.kg/s^2.K; Boltzmann constant
    g     = Pars.g;         % m/s^2; Gravitational constant
    H     = Pars.H;         % m; bed height
    Dint  = Pars.Dint;      % m; inside column diameter
    Dext  = Pars.Dext;      % m; outside column diameter
    rhow  = Pars.rhow;
    Cpw   = Pars.Cpw;
    kw    = Pars.kw;        % W/m.K; wall thermal conductivity
    ks    = Pars.ks;
    rhos  = Pars.rhos;
    Cps   = Pars.Cps;
    rhob  = Pars.rhob;
    dHads = Pars.dHads;
    Rp    = Pars.Dp;        % m; mean particle diameter
    rp    = Pars.rp;        % m; mean pore radius
    Tamb  = Pars.T0;        % K; ambient temperature
    T0    = Pars.T0;
    Tginlet = Pars.Tginlet;
    Patm  = Pars.Patm;      % Pa; atmospheric pressure (outside column)
    Ptot  = Pars.Ppa;
    cginlet = Pars.cginlet;
    vg    = Pars.vg;        % m/s; bulk interstitial velocity
    eb    = Pars.eb;        % --; bed void fraction = 0.375 + 0.34*Dp/Dint (Jeshar eqn)
    ep    = Pars.ep;        % --; particle void fraction
    Rint  = Dint/2;                   % m; radius of column 
    Rext  = Dext/2;                   % m; outer radius of column with a thickness of 9mm
    asint = Rint/(Rext^2 - Rint^2);
    asext = Rext/(Rext^2 - Rint^2);
    asp   = 3/Rp;
    %----------------------------
    tc = transferCoeffs(Tamb,Tamb,Pars);
    rhog = tc.rhog;
    Cpg = tc.Cpg;
    De = tc.De;
    Dz = tc.Dz;
    kf = tc.kf;
    Lz = tc.Lz;
    kz = tc.kz;
    hf = tc.hf;
    hw = tc.hi;
    Uo = tc.Uo;
    %----------------------------
    % Grid in axial direction
    nz  = 20;
    dz  = H/(nz - 1);
    dzs = dz^2;
    z   = [0:dz:H]';
    %z   = flip(log(linspace(exp(H),exp(0),nz)));
    
    % Grid in radial direction
    nr  = 20;
    if nr == 1
        dr = 0;
    else
        dr  = Rp/(nr - 1);
    end
    drs = dr^2;
    r   = [0:dr:Rp]';

    % Independant variable for ODE integration
    tf    = 120*60;                % s; final time
    nt    = 100;
    tout  = linspace(0,tf,nt)';

    % ------------------------------------------------------------
    % Initialize FD eqns
    cgz=zeros(nz,1); cgzz=zeros(nz,1); cgt=zeros(nz,1);
    Tgz=zeros(nz,1); Tgzz=zeros(nz,1); Tgt=zeros(nz,1);
    Twt=zeros(nz,1); 
    cpr=zeros(nz,nr); cprr=zeros(nz,nr); cpt=zeros(nz,nr);
    Tpr=zeros(nz,nr); Tprr=zeros(nz,nr); Tpt=zeros(nz,nr);
%     qct=zeros(nz,1); 

    cg    = zeros(nz,1); 
    cg(1) = cginlet;  
    Tg = zeros(nz,1)  + T0;  %Tg(:,1)=T0;
    Tw = zeros(nz,1)  + T0;  %Tw(:,1)=T0;
    cp = zeros(nz,nr);       %cp(:,:)=c0;
    Tp = zeros(nz,nr) + T0;  %Tp(:,:)=T0;
%     qc = zeros(nz,1);        %q(:,1)=T0;
    store = zeros(nz,nr);
    y0 = [cg; cp(:); Tg; Tp(:); Tw];

    % ODE integration
    reltol  = 1.0e-2; 
    abstol  = 1.0e-2;
    options = odeset('RelTol',reltol,'AbsTol',abstol);

    [t,y] = ode15s(@MOLScenerio1,tout,y0,options); %returns mat with nvar*nz*nr cols and t rows

    % Returns matrices that are ntxnz*nr size
    ysol.cg  = reshape(y(1:nt,1                : nz          ),nt,nz,1)/cginlet;
    ysol.cp  = reshape(y(1:nt,nz+1             : nz+nz*nr    ),nt,nz,nr);
    ysol.Tg  = reshape(y(1:nt,nz+nz*nr+1       : 2*nz+nz*nr  ),nt,nz,1);
    ysol.Tp  = reshape(y(1:nt,2*nz+nz*nr+1     : 2*nz+2*nz*nr),nt,nz,nr);
    ysol.Tw  = reshape(y(1:nt,2*nz+2*nz*nr+1   : 3*nz+2*nz*nr),nt,nz,1);
%     ysol.qc  = reshape(y(1:nt,3*nz+2*nz*nr+1 : 4*nz+2*nz*nr),nt,nz,1);
    
    ysol.t = t;
    ysol.y = y;
    ysol.z = z;
    ysol.r = r;
    ysol.dqdc = store;

    function yt = MOLScenerio1(t,y)
        % Step through the grid points in r and z:
        cg  = reshape(y(1              : nz          ),nz,1);
        cp  = reshape(y(nz+1           : nz+nz*nr    ),nz,nr);
        Tg  = reshape(y(nz+nz*nr+1     : 2*nz+nz*nr  ),nz,1);
        Tp  = reshape(y(2*nz+nz*nr+1   : 2*nz+2*nz*nr),nz,nr);
        Tw  = reshape(y(2*nz+2*nz*nr+1 : 3*nz+2*nz*nr),nz,1);
%         qc  = reshape(y(3*nz+2*nz*nr+1 : 4*nz+2*nz*nr),nz,1);
        
        for i=1:nz
            for j=1:nr       

            % ------------------------------------------------------------
            % FD eqns for r dependant variables
            if(j==1) 
                cp_bcw    = 0;
                Tp_bcw    = 0;
                cp_w      = cp(i,j+1) - (2*dz)*cp_bcw;
                Tp_w      = Tp(i,j+1) - (2*dz)*Tp_bcw;
                cpr(i,j)  = (cp(i,j+1) - 2*cp(i,j) + cp_w)/drs;     % cp(i,j-1)=cp(i,j)
                cprr(i,j) = (cp(i,j+1) - 2*cp(i,j) + cp_w)/drs;  
                Tpr(i,j)  = (Tp(i,j+1) - 2*Tp(i,j) + Tp_w)/drs;  
                Tprr(i,j) = (Tp(i,j+1) - 2*Tp(i,j) + Tp_w)/drs;  
            elseif(j==nr) 
                cp_bce    = (kf/De)*(cg(i)-cp(i,nr));
                Tp_bce    = (hf/ks)*(Tg(i)-Tp(i,nr));
                cp_e      = cp(i,j-1) + (2*dr)*cp_bce;
                Tp_e      = Tp(i,j-1) + (2*dr)*Tp_bce;
                cpr(i,j)  = (1/r(j))*(cp_e - cp(i,j-1))/(2*dr);     % cp(i,j-1)=cp(i,j)
                cprr(i,j) = (cp_e - 2*cp(i,j) + cp(i,j-1))/drs;  
                Tpr(i,j)  = (1/r(j))*(Tp_e - Tp(i,j-1))/(2*dr); 
                Tprr(i,j) = (Tp_e - 2*Tp(i,j) + Tp(i,j-1))/drs; 
            else 
                cpr(i,j)  = (1/r(j))*(cp(i,j+1) - cp(i,j-1))/(2*dr); 
                cprr(i,j) = (cp(i,j+1) - 2*cp(i,j) + cp(i,j-1))/drs;  
                Tpr(i,j)  = (1/r(j))*(Tp(i,j+1) - Tp(i,j-1))/(2*dr); 
                Tprr(i,j) = (Tp(i,j+1) - 2*Tp(i,j) + Tp(i,j-1))/drs; 
            end
            
            % ------------------------------------------------------------
            % FD eqns for z dependant variables
            if(i==1) 
                cg_bcw = -vg/Dz*(cginlet - cg(i));
                Tg_bcw = -rhog*Cpg*vg/Lz*(Tginlet - Tg(i));
                cg_w   = cg(i) - (dz)*cg_bcw;
                cg_e   = cg(i+1);
                Tg_w   = Tg(i) - (dz)*Tg_bcw;
                Tg_e   = Tg(i+1);
            elseif(i==nz) 
                cg_bce = 0;
                Tg_bce = 0;
                cg_w   = cg(i-1);
                cg_e   = cg(i) + (dz)*cg_bce;
                Tg_w   = Tg(i-1);
                Tg_e   = Tg(i) + (dz)*Tg_bce;
            else 
                cg_w = cg(i-1);
                cg_e = cg(i+1);
                Tg_w = Tg(i-1);
                Tg_e = Tg(i+1);
            end
            cgz(i)  = (cg(i) - cg_w)/dz;                    
            cgzz(i) = (cg_e - 2.0*cg(i) + cg_w)/dzs;    
            Tgz(i)  = (Tg(i) - Tg_w)/dz;                   
            Tgzz(i) = (Tg_e - 2.0*Tg(i) + Tg_w)/dzs;
            
            % ------------------------------------------------------------
            % PDEs
            % ------------------------------------------------------------                         
            %Toth = Toth_isotherm(cp(i,j),Tp(i,j));
            %IsoT     = DA_isotherm(cp(i,j),Tp(i,j),Ptot);
            IsoT   = LI_isotherm(cp(i,j),Tp(i,j));
            %IsoT   = Toth_isotherm(cp(i,j),Tp(i,j));
            dqpdt  = 3*kf/Rp/rhos*(cg(i) - cp(i,nr));
            dqdc   = IsoT.dqdc;
            %disp(dqdc);
            store(i,j) = dqdc;
            %dqdc   = 10;
            % DHads  = DA.DHads;
            alfa   = (1 + rhos*(1-ep)/ep*dqdc)^(-1);
            
            cgt(i)=...
                Dz                              *cgzz(i)...   Diffusive
                -vg                             *cgz(i)...    Convective
                -rhos*(1-eb)/eb*dqpdt...                                  Source
                ;

            Tgt(i)=...
                Lz/(eb*Cpg*rhog)                *Tgzz(i)...   Diffusive
                -vg                             *Tgz(i)...    Convective
                -asp*hf*(1-eb)/(eb*Cpg*rhog)    *(Tg(i)-Tp(i,nr))...
                -2*hw/(eb*Rint*Cpg*rhog)        *(Tg(i)-Tw(i))... Source
                ;               

            Twt(i)=...
                asint*hw/(rhow*Cpw)             *(Tg(i)-Tw(i))... Source
                -asext*Uo/(rhow*Cpw)            *(Tw(i)-Tamb)... Source
                ;
            
            cpt(i,j)=...
                alfa*De*2                       *cpr(i,j)...
                +alfa*De                        *cprr(i,j)...
                ;
            
            CPe = 1/(ep*rhog*Cpg + (1 - ep)*rhos*Cps);
            Tpt(i,j)=...
                CPe*ks                          *Tprr(i,j)...
                +CPe*ks*2                       *Tpr(i,j)...
                +CPe*(dHads)*rhob*cpt(i,j)*dqdc...                 
                ;
            
            end
        end
        % ------------------------------------------------------------
        % 2D to 1D matrices
        yt    = [cgt; cpt(:); Tgt; Tpt(:); Twt];

    end
    
end

function LI = LI_isotherm(c,T)
    Rg = 8.314;             % J/mol.K
    Mv = 18.01528e-3;       % kg/mol
    
    P  = c*Rg*T;
    k1 = 174.2;             % mol/kg
    k2 = 0.510;             % mol/kg.K
    k3 = 1.375/101325;      % 1/Pa; 1/atm -> 1/Pa (1/101325)
    %k4 = 408.2;             % K
    k4 = 310.2;  
    
    %qm = (k1 - k2*T)*Mv;    % mol/kg -> kg_v/kg_s
    qm = 30*(k1 - k2*T);    % mol/kg -> kg_v/kg_s
    b  = k3*exp(k4/T);      % 1/Pa
    
    dqdP     = (b*qm)/(b*P + 1)^2;
    LI.dqdc  = Rg*T*dqdP;
end

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

    P = cp*R*Tp;
    Toth.qce = a*P/(1 + (b*P)^n)^(1/n);
    dqdP = a*((b*P)^n + 1)^(-(n+1)/n);
    Toth.dqdc = 1e-3*R*Tp*dqdP/2;
end

function DA = DA_isotherm(c,T,Ptot)
    Rg = 8.314; % J/mol.K
    Mv = 18.01528e-3; % kg/mol
    
    Pv    = c*Rg*T;
    Psv   = CoolProp.PropsSI('P','T',T,'Q',0,'Water');
    lnPr  = log(Psv/Pv);
    
    Hsv = CoolProp.PropsSI('H','T',T,'Q',1,'Water');
    Hsl = CoolProp.PropsSI('H','T',T,'Q',0,'Water');
    
    DHvap = Hsl - Hsv;
    % DHvap1 = -7.1*T^2 + 2620.4*T + 2.2866e6; % J/kg
    
    beta = CoolProp.PropsSI('isobaric_expansion_coefficient','P',Ptot,'T',T,'Water');
    beta20 = CoolProp.PropsSI('isobaric_expansion_coefficient','P',Ptot,'T',273.15+20,'Water');
    rho20  = CoolProp.PropsSI('Dmass','P',Ptot,'T',273.15+20,'Water');
    rho_ads = rho20/(1 + beta20*(T - 293.15));
    
    Wo  = 341.03;                % ml/g -> m3/kg;
    E   = 1192.3e3;              % J/kg;
    E2  = E*Mv;                  % J/mol;
    n   = 1.55;                  % --;
    RTE = Rg*T/E2;
    
    A  = Rg/Mv*T*lnPr; 
    
    W  = Wo*exp(-(A/E)^n);      % m3/kg
    theta = W/Wo;
    qs = Wo*rho_ads/Mv;         % mol/kg
    q  = W*rho_ads/Mv;          % mol/kg
    
    dqdP_T = n*qs*exp(-RTE^n*lnPr^n)*RTE^n*lnPr^n/(Pv*lnPr);
    dqdc   = Rg*T*dqdP_T;       % (mol_v/kg_s)/(mol_v/m3_v) = m3_v/kg_s
    
    DH   = DHvap + A;
    DH2  = DHvap + E*(log(1/theta))^(1/n) + E*beta*T/n*(log(1/theta))^-(n-1)/n; 

    DA.dqdc = dqdc;
    DA.DH   = DH;
    DA.DH2  = DH2;
end


