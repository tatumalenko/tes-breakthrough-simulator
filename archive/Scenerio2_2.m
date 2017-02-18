function ysol = Scenerio2_2(varargin)
    % ---------------------------------------------------------------------
    % PROCESS VARIABLE INPUT ARGUMENTS
    % ---------------------------------------------------------------------
    Pars = Pars_TES();
    
    isoA = 1;
    isoN = 1;
    if nargin > 1
        for iLo = 1:2:nargin-1
            if strcmp(varargin{iLo},'isotherm');
                isotherm = varargin{iLo+1};
            elseif strcmp(varargin{iLo},'dH')
                dHads = varargin{iLo+1};   
            elseif strcmp(varargin{iLo},'isoA')
                isoA = varargin{iLo+1};  
            elseif strcmp(varargin{iLo},'isoN')
                isoN = varargin{iLo+1};   
            end
        end
    end
    
    if isempty(varargin)
        isotherm = 30;
    end
    
    % ---------------------------------------------------------------------
    % IMPORT MODEL PARAMETERS (CONSTANTS)
    % ---------------------------------------------------------------------
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
    rhop  = Pars.rhop;
    %dHads = Pars.dHads;
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
    et    = eb + (1 - eb)*ep;


    % ---------------------------------------------------------------------
    % GRID IN AXIAL DIRECTION
    % ---------------------------------------------------------------------
    nz  = 20;
    dz  = H/(nz - 1);
    dzs = dz^2;
    z   = [0:dz:H]';
    
    % ---------------------------------------------------------------------
    % GRID IN RADIAL DIRECTION
    % ---------------------------------------------------------------------
    nr  = 20;
    if nr == 1
        dr = 0;
    else
        dr  = Rp/(nr - 1);
    end
    drs = dr^2;
    r   = [0:dr:Rp]';

    % ---------------------------------------------------------------------
    % INDEPENDENT VARIABLE FOR ODE INTEGRATION
    % ---------------------------------------------------------------------
    tf    = 120*60;                % s; final time
    nt    = 100;
    tout  = linspace(0,tf,nt)';

    % ---------------------------------------------------------------------
    % INITIALIZE VARIABLES
    % ---------------------------------------------------------------------
    cgz=zeros(nz,1); cgzz=zeros(nz,1); cgt=zeros(nz,1);
    Tgz=zeros(nz,1); Tgzz=zeros(nz,1); Tgt=zeros(nz,1);
    Twt=zeros(nz,1); 
    cpr=zeros(nz,nr); cprr=zeros(nz,nr); cpt=zeros(nz,nr);
    Tpr=zeros(nz,nr); Tprr=zeros(nz,nr); Tpt=zeros(nz,nr);
    %qct=zeros(nz,1); 
    HADS = zeros(nz,nr);
    
    % ---------------------------------------------------------------------
    % SET INITIAL CONDITIONS
    % ---------------------------------------------------------------------
    cg = zeros(nz,1); 
    Tg = zeros(nz,1)  + T0;  %Tg(:,1)=T0;
    Tw = zeros(nz,1)  + T0;  %Tw(:,1)=T0;
    cp = zeros(nz,nr);       %cp(:,:)=c0;
    Tp = zeros(nz,nr) + T0;  %Tp(:,:)=T0;
    %qc = zeros(nz,1);        %q(:,1)=T0;

    y0 = [cg; cp(:); Tg; Tp(:); Tw];

    % ---------------------------------------------------------------------
    % ODE INTEGRATION OPTIONS
    % ---------------------------------------------------------------------
    reltol  = 1.0e-2; 
    abstol  = 1.0e-2;
    options = odeset('RelTol',reltol,'AbsTol',abstol);

    [t,y] = ode15s(@MOLScenerio1,tout,y0,options); %returns mat with nvar*nz*nr cols and t rows

    % ---------------------------------------------------------------------
    % RETURNS MATRICES THAT ARE ntxnz*nr IN SIZE
    % ---------------------------------------------------------------------
    ysol.cg  = reshape(y(1:nt,1                : nz          ),nt,nz,1);
    ysol.cp  = reshape(y(1:nt,nz+1             : nz+nz*nr    ),nt,nz,nr);
    ysol.Tg  = reshape(y(1:nt,nz+nz*nr+1       : 2*nz+nz*nr  ),nt,nz,1);
    ysol.Tp  = reshape(y(1:nt,2*nz+nz*nr+1     : 2*nz+2*nz*nr),nt,nz,nr);
    ysol.Tw  = reshape(y(1:nt,2*nz+2*nz*nr+1   : 3*nz+2*nz*nr),nt,nz,1);
    %ysol.qc  = reshape(y(1:nt,3*nz+2*nz*nr+1 : 4*nz+2*nz*nr),nt,nz,1);
    
    ysol.t       = t;
    ysol.y       = y;
    ysol.z       = z;
    ysol.r       = r;
    ysol.cginlet = cginlet;
    ysol.dH      = HADS;
    ysol.cgend   = ysol.cg(:,end)./cginlet;
    ysol.isotherm= isotherm;
    ysol.isoA    = isoA;
    ysol.isoN    = isoN;
    ysol.dqdc    = isotherm.*(1-(isoA.*(ysol.cgend)).^isoN);

    function yt = MOLScenerio1(t,y)
        % -----------------------------------------------------------------
        % UNPACK Y (A SINGLE ROW VECTOR) INTO EACH DEPENDENT VARIABLE
        % -----------------------------------------------------------------
        cg  = reshape(y(1              : nz          ),nz,1);
        cp  = reshape(y(nz+1           : nz+nz*nr    ),nz,nr);
        Tg  = reshape(y(nz+nz*nr+1     : 2*nz+nz*nr  ),nz,1);
        Tp  = reshape(y(2*nz+nz*nr+1   : 2*nz+2*nz*nr),nz,nr);
        Tw  = reshape(y(2*nz+2*nz*nr+1 : 3*nz+2*nz*nr),nz,1);
        %qc  = reshape(y(3*nz+2*nz*nr+1 : 4*nz+2*nz*nr),nz,1);

        cg(cg<0) = 0.0;
        cp(cp<0) = 0.0;
        
        assert(all(isreal(cg(:)))&&all(isreal(cp(:)))&&all(cg(:)>=0)&&all(cp(:)>=0));
        
        % -----------------------------------------------------------------
        % STEP THROUGH THE GRID POINTS IN r AND z:
        % -----------------------------------------------------------------
        for i=1:nz
            % -----------------------------------------------------------------
            % IMPORT TRANSFER COEFFICIENTS (TEMP DEPENDENT)
            % -----------------------------------------------------------------
            tc   = transferCoeffs(Tg(i),Tw(i),Pars);
            rhog = tc.rhog;
            Cpg  = tc.Cpg;
            De   = tc.De;
            Dz   = tc.Dz;
            kf   = tc.kf;
            Lz   = tc.Lz;
            kz   = tc.kz;
            hf   = tc.hf;
            hw   = tc.hi;
            ho   = tc.ho;
            Uo   = tc.Uo;
            kg   = tc.kg;
            
            for j=1:nr       
            % -------------------------------------------------------------
            % FD eqns for r dependant variables
            % -------------------------------------------------------------
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
            
            % -------------------------------------------------------------
            % FD eqns for z dependant variables
            % -------------------------------------------------------------
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
            
            % -------------------------------------------------------------
            % ISOTHERM SELECTION
            % -------------------------------------------------------------
            if ~isnumeric(isotherm)
                if strcmp(isotherm,'DA')
                    % IsoT   = DA_isotherm(cp(i,j),Tp(i,j),Ptot);
                    % dqdc   = IsoT.dqdc;
                    % DHads  = IsoT.DHads;
                elseif strcmp(isotherm,'LI')
                    IsoT   = LI_isotherm(cp(i,j),Tp(i,j));
                    dqdc   = IsoT.dqdc;
                    dHads  = IsoT.dh;
                elseif strcmp(isotherm,'Toth')
%                     Iso1   = Toth_isotherm(cp(i,j),Tp(i,j));
%                     Iso2   = Toth_isotherm(cp(i,j-1),Tp(i,j));
%                     dqdc   = (Iso1.q - Iso2.q)/(cp(i,j) - cp(i,j-1));
                    IsoT   = Toth_isotherm(cp(i,j),Tp(i,j));
                    dqdc   = IsoT.dqdc;
                    dHads  = 82000;
                elseif strcmp(isotherm,'LI')
%                     IsoT   = LI_isotherm(cp(i,j),Tp(i,j));
%                     dqdc   = IsoT.dqdc;
                    dHads  = 75000;
                    Rgg  = R_GAS('kPa','m3');   % kPa.m3/molG.K
                    P  = c*Rgg*T;             % molG/m3G*kPa.m3G/molG/K*K = kPa

                    k1 = 24.88;            % molA/kgS
                    k2 = 0.0422;           % molA/kgS.K
                    k3 = 1.3550e-05;       % 1/kPa
                    k4 = 4865;             % K 

                    qm = (k1 - k2*T);      % molA/kgS
                    b  = k3*exp(k4/T);   % 1/kPa
                    dqdp     = (qm*b)/(1 + b*P)^2;   % molA/kgS.kPa
                    dqdc     = dqdp*Rgg*T;               % (molA/molG)*m3G/kgS = molA/kgS*m3G/molG 

                    if dqdc > 100
                        dqdc = 100;
                    elseif dqdc < 1
                        dqdc = 1;
                    end
                end              
            else
                dqdc   = isotherm;
%                 dHads  = IsoT.dh;
%                 norm   = cg(i)./cginlet;
%                 dqdc   = isotherm*(1 - (isoA*norm)^isoN);
                dHads  = 75000;
            end
       
            HADS(i,j) = dHads;
%             if dqdc < 1 || dqdc > 150
%                 disp(dqdc);
%             end

            % -------------------------------------------------------------
            % PDES
            % -------------------------------------------------------------
            dqpdt  = 3*kf/Rp/rhos*(cg(i) - cp(i,nr));
            alfa   = (1 + rhos*(1-ep)/ep*dqdc)^(-1);
            
            cgt(i)=...
                Dz                              *cgzz(i)...                 Diffusive
                -vg                             *cgz(i)...                  Convective
                -rhos*(1-eb)/eb*dqpdt...                                    Source
                ;

            Tgt(i)=...
                kg/(eb*Cpg*rhog)                *Tgzz(i)...                 Diffusive
                -vg                             *Tgz(i)...                  Convective
                -3/Rp*hf*(1-eb)/(eb*Cpg*rhog)   *(Tg(i)-Tp(i,nr))...
                -2/Rint*hw/(eb*Cpg*rhog)        *(Tg(i)-Tw(i))...           Source
                ;               

            Twt(i)=...
                2*Rint*hw/((Rext^2 - Rint^2)*rhow*Cpw)     *(Tg(i)-Tw(i))...           Source
                -2*Rext*ho/((Rext^2 - Rint^2)*rhow*Cpw)    *(Tw(i)-Tamb)...            Source
                ;
            
            cpt(i,j)=...
                alfa*De*2                       *cpr(i,j)...
                +alfa*De                        *cprr(i,j)...
                ;
            
            CPe = 1/(ep*rhog*Cpg + (1 - ep)*rhos*Cps);
            Tpt(i,j)=...
                CPe*ks                          *Tprr(i,j)...
                +CPe*ks*2                       *Tpr(i,j)...
                +CPe*(dHads)*rhop*cpt(i,j)*dqdc...                 
                ;
            
            end
        end
        % -----------------------------------------------------------------
        % 2D to 1D matrices
        % -----------------------------------------------------------------
        yt    = [cgt; cpt(:); Tgt; Tpt(:); Twt];

    end
    
end








