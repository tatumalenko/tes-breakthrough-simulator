function ysol = breakthrough_scenario
% ---------------------------------------------------------------------
% PROCESS VARIABLE INPUT ARGUMENTS
% ---------------------------------------------------------------------

% model_option
% 1: instantaneous adsorption 1-D
% cl----> cl_t = Dz*cl_zz - ul*cl_z - src/ec
% qads--> qads_t = src/((1-ec)*rhos)
% 2: instantaneous adsorption 2-D
% cl----> cl_t = Dz*cl_zz - ul*cl_z - src_bed/ec
% qads--> cp_t = ep*De*(2/r*cp_r + cp_rr)/(ep + (1-ep)*rhos*dqdcp)
% 3: adsorption kinetic limiting 2-D
% cl----> cl_t = Dz*cl_zz - ul*cl_z - src_bed/ec
% cp----> cp_t = De*(2/r*cp_r + cp_rr) - src_pel/ep
% qads--> qads = src_pel/((1 - ep)*rhos)

% adsorption_option
% 1: fluid phase DF (clstar,cpstar)
% src = kl*abed*(cl - clstar)
% src_bed = kl*abed*(cl - cpend)
% src_pel = kads*apel*(cp - cpstar)
% 2: solid phase DF (castar)
% src = kl*abed*(castar - ca)
% src_bed = kl*abed*(cl - cpend)
% src_pel = kads*apel*(cp - cpstar)

% isotherm_option
% 1: 


% ---------------------------------------------------------------------
% IMPORT MODEL PARAMETERS (CONSTANTS)
% ---------------------------------------------------------------------
PI            = pi();
[~,~,simpars] = xlsread('/Users/patriceamyot/Desktop/Breakthrough 20 Oct 2015 - Patrice - REV03.xlsm','Sheet1','B2:B14');
[~,~,Sheet1]  = xlsread('/Users/patriceamyot/Desktop/Breakthrough 20 Oct 2015 - Patrice - REV03.xlsm','Sheet1');
[~,~,Sheet2]  = xlsread('/Users/patriceamyot/Desktop/Breakthrough 20 Oct 2015 - Patrice - REV03.xlsm','Sheet2');
expt.t        = cell2mat(Sheet2(4:21,1));
expt.cl       = cell2mat(Sheet2(4:21,2:8));

model_option      = simpars{1};
tf                = simpars{2};
nz                = 30; %simpars{3};
nr                = simpars{4};
flow              = simpars{5};
Hbed              = simpars{6};
Rbed              = simpars{7};
eb                = simpars{8};
Rpel              = simpars{9};
mass              = simpars{10};
rhos              = simpars{11};
rowiso            = simpars{12};
adsorption_option = simpars{13};

if rowiso == 1
    irow = 20;
    nc   = Sheet1{IROW,2};
    isotherm_option = Sheet1{IROW,5};
    % Read information on adsorbent && isotherms
    for(i = 1:nc)
        namecomp{i} = Sheet1{irow + 1 + i,1};
        clinlet(i)  = Sheet1{irow + 1 + i,3};
        pads(i,1)   = Sheet1{irow + 1 + i,4}/1000;  % Changed from g/kg:kg/kg
        pads(i,2)   = Sheet1{irow + 1 + i,5};
        pads(i,3)   = Sheet1{irow + 1 + i,6};
        pads(i,4)   = Sheet1{irow + 1 + i,7};
        kl(i)       = Sheet1{irow + 1 + i,8};
        Dz(i)       = Sheet1{irow + 1 + i,9};
        De(i)       = Sheet1{irow + 1 + i,10};
        kads(i)     = Sheet1{irow + 1 + i,11};
    end
elseif rowiso == 2
    irow = 31;
    nc   = Sheet1{irow,2};
    isotherm_option = Sheet1{irow,5};
    % Read information on adsorbent && isotherms
    for i = 1:nc
        namecomp{i} = Sheet1{irow + 1 + i,1};
        clinlet(i)  = Sheet1{irow + 1 + i,3};
        for j = 1:9
            pads(i,j) = Sheet1{irow + 1 + i,3 + j};
        end
        kl(i)   = Sheet1{irow + 1 + i,13};
        Dz(i)   = Sheet1{irow + 1 + i,14};
        De(i)   = Sheet1{irow + 1 + i,15};
        kads(i) = Sheet1{irow + 1 + i,16};
    end
    pads(:,1) = pads(:,1)./1000;
else
    return;
end

% Some calculations    
rhob                = mass/(PI*Rbed^2*Hbed); % [kgS/m3B]
rhop                = rhob/(1 - eb);         % [kgS/m3S = kgS/m3B.m3B/m3S]
ep                  = 1 - rhop/rhos;
ulsup               = flow/(PI*Rbed^2);
ul                  = ulsup/eb;
abed                = (1 - eb)*(3/Rpel);
apel                = 250;

% ---------------------------------------------------------------------
% GRID IN AXIAL DIRECTION
% ---------------------------------------------------------------------
% nz  = simpars{3};
dz  = Hbed/(nz - 1);
dzs = dz^2;
z   = [0:dz:Hbed]';

% ---------------------------------------------------------------------
% GRID IN RADIAL DIRECTION
% ---------------------------------------------------------------------
if(model_option == 1)
    nr  = 1;
else
    nr  = 20; %simpars{4};
end

if nr == 1
    dr = 0;
else
    dr  = Rpel/(nr - 1);
end
drs = dr^2;
r   = [0:dr:Rpel]';

% ---------------------------------------------------------------------
% INDEPENDENT VARIABLE FOR ODE INTEGRATION
% ---------------------------------------------------------------------
% tf   = simpars{2};                % s; final time
nt    = 200;
tout  = linspace(0,tf,nt)';

% ---------------------------------------------------------------------
% INITIALIZE VARIABLES
% ---------------------------------------------------------------------
clz=zeros(nc,nz); clzz=zeros(nc,nz); clt=zeros(nc,nz);
cpr=zeros(nc,nz,nr); cprr=zeros(nc,nz,nr); cpt=zeros(nc,nz,nr);
cat=zeros(nc,nz,nr);

% ---------------------------------------------------------------------
% SET INITIAL CONDITIONS
% ---------------------------------------------------------------------
cl = zeros(nc,nz);
cp = zeros(nc,nz,nr);       %cp(:,:)=c0;
ca = zeros(nc,nz,nr);       %ca(:,:)=T0;

if(model_option == 1)
    y0 = [cl(:); ca(:)];
elseif(model_option == 2)
    y0 = [cl(:); cp(:)];
elseif(model_option == 3)
    y0 = [cl(:); cp(:); ca(:)];
else
    error('Invalid model_option number');
end

% ---------------------------------------------------------------------
% ODE INTEGRATION OPTIONS
% ---------------------------------------------------------------------
reltol  = 1.0e-2;
abstol  = 1.0e-2;
options = odeset('RelTol',reltol,'AbsTol',abstol,'NonNegative',[1:length(y0)]);

[t,y] = ode15s(@model,tout,y0,options); %returns mat with nvar*nz*nr cols and t rows

% ---------------------------------------------------------------------
% RETURNS MATRICES THAT ARE ntxnz*nr IN SIZE
% ---------------------------------------------------------------------
if(model_option == 1) 
    ysol.cl  = reshape(y(1:nt,1       : nc*nz          ),nt,nc,nz);
    ysol.ca  = reshape(y(1:nt,nc*nz+1 : nc*nz+nc*nz*nr ),nt,nc,nz,nr);
elseif(model_option == 2)
    ysol.cl  = reshape(y(1:nt,1       : nc*nz          ),nt,nc,nz);
    ysol.cp  = reshape(y(1:nt,nc*nz+1 : nc*nz+nc*nz*nr ),nt,nc,nz,nr);
elseif(model_option == 3)
    ysol.cl  = reshape(y(1:nt,1                : nc*nz             ),nt,nc,nz);
    ysol.cp  = reshape(y(1:nt,nc*nz+1          : nc*nz+nc*nz*nr    ),nt,nc,nz,nr);
    ysol.ca  = reshape(y(1:nt,nc*nz+nc*nz*nr+1 : nc*nz+2*nc*nz*nr  ),nt,nc,nz,nr);
end

ysol.t       = t;
ysol.y       = y;
ysol.z       = z;
ysol.r       = r;
ysol.clinlet = clinlet;
ysol.clexit  = ysol.cl(:,:,end);

    function yt = model(t,y)
        % -----------------------------------------------------------------
        % UNPACK Y (A SINGLE ROW VECTOR) INTO EACH DEPENDENT VARIABLE
        % -----------------------------------------------------------------
        if(model_option == 1)
            cl  = reshape(y(1                 : nc*nz            ),nc,nz);
            ca  = reshape(y(nc*nz+1           : nc*nz+nc*nz*nr   ),nc,nz,nr);
        elseif(model_option == 2)
            cl  = reshape(y(1                 : nc*nz            ),nc,nz);
            cp  = reshape(y(nc*nz+1           : nc*nz+nc*nz*nr   ),nc,nz,nr);
        elseif(model_option == 3)
            cl  = reshape(y(1                 : nc*nz            ),nc,nz);
            cp  = reshape(y(nc*nz+1           : nc*nz+nc*nz*nr   ),nc,nz,nr);
            ca  = reshape(y(nc*nz+nc*nz*nr+1  : nc*nz+2*nc*nz*nr ),nc,nz,nr);
        end
        
%         cl(:,1) = clinlet(:);
        cl(cl<0) = 0.0; 
        cp(cp<0) = 0.0;
        ca(ca<0) = 0.0;
        
        
%         assert(all(isreal(cl(:)))&&all(isreal(cp(:)))&&all(isreal(ca(:)))...
%             &&all(cl(:)>=0)&&all(cp(:)>=0)&&all(ca(:)>=0));
        
        % -----------------------------------------------------------------
        % STEP THROUGH THE GRID POINTS IN r AND z:
        % -----------------------------------------------------------------
        for ic = 1:nc
            for iz = 1:nz
                % -----------------------------------------------------------------
                % IMPORT TRANSFER COEFFICIENTS (TEMP DEPENDENT)
                % -----------------------------------------------------------------

                for ir = 1:nr
                    % -------------------------------------------------------------
                    % FD eqns for r dependant variables
                    % -------------------------------------------------------------
                    if(nr > 1)
                        if(ir == 1)
                            cp_bcw         = 0;
                            cp_w           = cp(ic,iz,ir+1) - (2*dz)*cp_bcw;
                            cpr(ic,iz,ir)  = (cp(ic,iz,ir+1) - 2*cp(ic,iz,ir) + cp_w)/drs;     % cp(i,j-1)=cp(i,j)
                            cprr(ic,iz,ir) = (cp(ic,iz,ir+1) - 2*cp(ic,iz,ir) + cp_w)/drs;
                        elseif(ir == nr)
                            cp_bce         = (kl(ic)/De(ic))*(cl(ic,iz)-cp(ic,iz,nr));
                            cp_e           = cp(ic,iz,ir-1) + (2*dr)*cp_bce;
                            cpr(ic,iz,ir)  = (1/r(ir))*(cp_e - cp(ic,iz,ir-1))/(2*dr);     % cp(i,ir-1)=cp(i,ir)
                            cprr(ic,iz,ir) = (cp_e - 2*cp(ic,iz,ir) + cp(ic,iz,ir-1))/drs;
                        else
                            cpr(ic,iz,ir)  = (1/r(ir))*(cp(ic,iz,ir+1) - cp(ic,iz,ir-1))/(2*dr);
                            cprr(ic,iz,ir) = (cp(ic,iz,ir+1) - 2*cp(ic,iz,ir) + cp(ic,iz,ir-1))/drs;
                        end
                    end
                    % -------------------------------------------------------------
                    % FD eqns for z dependant variables
                    % -------------------------------------------------------------
                    if(iz == 1)
%                         cl_bcw = 0;
                        cl_bcw = -ul/Dz(ic)*(clinlet(ic) - cl(ic,iz));
                        cl_w   = cl(ic,iz) - (dz)*cl_bcw;
                        cl_e   = cl(ic,iz+1);
                    elseif(iz == nz)
                        cl_bce = 0;
                        cl_w   = cl(ic,iz-1);
                        cl_e   = cl(ic,iz) + (dz)*cl_bce;
                    else
                        cl_w = cl(ic,iz-1);
                        cl_e = cl(ic,iz+1);
                    end
                    clz(ic,iz)  = (cl(ic,iz) - cl_w)/dz;
                    clzz(ic,iz) = (cl_e - 2.0*cl(ic,iz) + cl_w)/dzs;
                    
                    % -------------------------------------------------------------
                    % PDES
                    % -------------------------------------------------------------
                    
                    switch model_option
                        case 1 % 1-D LDF Model (either with film DF or solid DF)
                            if adsorption_option == 1
                                clstar = calcclstar(ic,iz,ir,ca,pads,isotherm_option);
                                src    = kl(ic)*(3/Rpel)*(cl(ic,iz) - clstar);
                            elseif adsorption_option == 2
                                castar = calccastar(ic,iz,ir,cl,pads,isotherm_option);
                                keff   = 15*De(1)/Rpel^2;
                                src    = keff*(castar - ca(ic,iz,ir));
                            end
                            
                            clt(ic,iz)=...
                                Dz(ic)                *clzz(ic,iz)...   Diffusive
                                -ul                   *clz(ic,iz)...    Convective
                                -src*rhob*(1-eb)/eb...                  Source
                                ;
                            
                            cat(ic,iz,ir)=...
                                src;
                           
                        case 2
                            src_bed = kl(ic)*abed*(cl(ic,iz) - cp(ic,iz,end));
                            dqdcp   = calcdqdcp(cp(ic,iz,ir),isotherm_option);
                            alfa    = (ep + (1-ep)*rhos*dqdcp)^(-1);
                            
                            clt(ic,iz)=...
                                Dz(ic)                *clzz(ic,iz)...   Diffusive
                                -ul                   *clz(ic,iz)...    Convective
                                -src_bed/eb...                          Source
                                ;
                            
                            cpt(ic,iz,ir)=...
                                ep*De(ic)*alfa*2      *cpr(ic,iz,ir)...
                                +ep*De(ic)*alfa       *cprr(ic,iz,ir)...
                                ;
                         
                        case 3
                            if adsorption_option == 1
                                cpstar = calccpstar(ic,iz,ir,ca,pads,isotherm_option);
                                src_pel = kads(ic)*apel*(cp(ic,iz,ir) - cpstar);
                            elseif adsorption_option == 2
                                castar = calccastar(ic,iz,ir,cp,pads,isotherm_option);
                                src_pel = kads(ic)*apel*(castar - ca(ic,iz,ir));
                            end
                            
                            src_bed = kl(ic)*abed*(cl(ic,iz) - cp(ic,iz,end));
                            
                            clt(ic,iz)=...
                                Dz(ic)                *clzz(ic,iz)...   Diffusive
                                -ul                   *clz(ic,iz)...    Convective
                                -src_bed/eb...                          Source
                                ;
                            
                            cpt(ic,iz,ir)=...
                                De(ic)*2             *cpr(ic,iz,ir)...
                                +De(ic)              *cprr(ic,iz,ir)...
                                -src_pel/ep...
                                ;
                            
                            cat(ic,iz,ir)=...
                                src_pel/((1-ep)*rhos);
                    end
                end
            end
        end
        if model_option == 1
            yt = [clt(:); cat(:)];
        elseif model_option == 2
            yt = [clt(:); cpt(:)];
        elseif model_option == 3
            yt = [clt(:); cpt(:); cat(:)];
        end
    end

end

function castar = calccastar(ic,iz,ir,conc,pads,isotherm_option)
nc = size(conc,1);
if length(size(conc)) == 2 % cl was sent in as conc
    if iz ~= 1
        concold = conc(:,iz-1);
    else
        concold = 0;
    end
    conc = conc(:,iz);
else                       % cp was sent in as conc
    conc = conc(:,iz,ir);
end

switch isotherm_option
    case 1 % Sips
        return
        
    case 2 % Freundlich
        % qi = KiCi(ai1C1 + ai2C2 + ai3C3 + ...)^(1/n1 - 1)
        AIJCJ    = 0;
        AIJCJold = 0;
        for ij = 1:nc          
            AIJCJ    = AIJCJ + pads(ic,ij+2)*conc(ij);
            if iz ~= 1
                AIJCJold = AIJCJold + pads(ic,ij+2)*concold(ij);
            end
        end
        if AIJCJ ~= 0
            castar = pads(ic,1)*conc(ic)*(AIJCJ)^(1/pads(ic,2) - 1);
%         elseif AIJCJold ~= 0
%             warning('calccastar: AIJCJold used in calculation');
%             castar = pads(ic,1)*concold(ic)*(AIJCJold)^(1/pads(ic,2) - 1);
        else
            castar = 0;
            %error('calccastar: Division by zero when all(conc)=0');
        end
        
    case 3 %
        return
end
end




