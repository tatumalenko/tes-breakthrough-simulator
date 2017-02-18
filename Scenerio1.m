%% *1. Scenerio1 function definition*
function ysol = Scenerio1(varargin)
    
    global T0 cginlet tao derr PI vi Dz % declare global variables
    PI = pi();
    
    nargs = length(varargin); % counts number arguments passed in
    iXLS = 0; % initialize the position of xls filename in argument to zero
    iWS = 0; % " " " of xls sheet name in argument to zero, WS=worksheet shorthand
    
    if isempty(varargin) % if empty, then no xls output and set default mesh sizes
        tM = 100;
        xN = 20;
    else                 % else, check position of name-value pair args, if
                         % any, and apply user-supplied mesh sizes, if any.
        for i=1:nargs
            if strcmp(varargin{i},'xlsname')
                xlsFile = varargin{i+1};
                iXLS = i;
            end

            if strcmp(varargin{i},'sheet')
                sheetName = varargin{i+1};
                iWS = i;
            end

            if iXLS==1||iWS==1
                tM = 1000;
                xN = 100;
            else
                tM = varargin{1};
                xN = varargin{2};
            end
        end
    end

    c_gExp24 = readtable('Dansdata.xlsx','Sheet','24 LPM Conc',...
                        'ReadVariableName',false);
    T_gExp24 = readtable('Dansdata.xlsx','Sheet','24 LPM Temp',...
                        'ReadVariableName',false);
    
    derr = 1e-15;   
    m = 0; % defines geometry of space mesh (0=rectangular)
    H = 0.0695; % m; height of column 
    t_f = 8000; % s; final time
    x = linspace(0,H,xN); % creates a linearly spaced array from 0 to H 
                          % divided in xN equally sized segments
    t = linspace(0,t_f,tM); % same idea as above but for a time array

    % Call pdepe MATLAB function pdepe with functions created (pdefun,
    % pdeic, pdebc) along with both x and t arrays defining both meshes and
    % assign the soln to variable 'sol'. NOTE: sol (the output of pdepe) is
    % an (tM by xN by K) array (i.e. a 3D array), where tM=no. points in
    % time mesh, xN=no. points in space mesh, K=no. of dependant variables
    % characterizing each PDEs as defined in the output variables of our
    % created function pdefun, where
    % u=4x1array=[u(1),u(2),u(3),u(4)]=[cg,Tg,Tw,q]-->see pdefun!
    try
        sol = pdepe(m,@pdefun,@pdeic,@pdebc,x,t);
    catch ME
        sol = pdepe(m,@pdefun,@pdeic,@pdebc,x,t);
    end
    
   
    c_g = sol(:,:,1); % the 2D array for conc of water in gas in mol/m3
    T_g = sol(:,:,2)-273.15; % the 2D array for bulk gas temp in C
    T_w = sol(:,:,3)-273.15; % the 2D array for inner wall temp in C
    q_p = sol(:,:,4); % the 2D array for solid's water capacity in kg/kg
    
    ysol.t  = t;
    ysol.cg = c_g; 
    ysol.Tg = T_g;
    ysol.Tw = T_w;
    ysol.q  = q_p;
  
    tCell = cell(tM,1);
    cg_outlet = cell(tM,1);
    Tg_outlet = cell(tM,1);
    Tw_outlet = cell(tM,1);
    qp_outlet = cell(tM,1);
    
   
    for i=1:tM
        tCell{i} = t(i);
        cg_outlet{i} = c_g(i,end); % breakthrough curve @ x=H (end=last index)
        Tg_outlet{i} = T_g(i,end); 
        Tw_outlet{i} = T_w(i,end);
        qp_outlet{i} = q_p(i,end);
    end


    colNames = {'Time, s','Cg(t,x=H), mol/m3','Tg(t,x=H), C',...
                'Tw(t,x=H), C','q(t,x=H), kg/kg'};

  
    colVars = [tCell,cg_outlet,Tg_outlet,Tw_outlet,qp_outlet];

   
    xlsMatrix = [colNames; colVars];
    
    if iXLS ~= 0 % if a xlsfile is specified then...
        if iWS == 0 % default first sheet of xlsfile
           xlswrite(xlsFile,xlsMatrix);
        else % else writes to specified sheet name
           xlswrite(xlsFile,xlsMatrix,sheetName); 
        end
    end
    
    createGUI;
    figure; 
    subplot(1,2,1); plot(c_gExp24.Var1,c_gExp24.Var2,'o',t,c_g(:,end)./cginlet); ylim([0 1]);
    subplot(1,2,2); plot(T_gExp24.Var1,T_gExp24.Var2-273.15,'o',t,T_g(:,end)); ylim([20 80]);
    %% *1A. pdefun function definition*
    function [c,f,s] = pdefun(x,t,u,DuDx)
        % Note to myself: Try fix, if numel(u)==2, call Scenerio1 again.
        R1 = 0.008313; % m3.kPa/mol.K
        R = 0.082; % L.atm/mol.K; ideal gas law constant
        T0 = 298.15; % K; initial/surrounding/inlet temperature
        Tav = (u(3)+T0)/2;
        Patm = 1.1; % atm; operating pressure
        Ppa = Patm*101325; % Pa; operating Pressure 
        yw = 0.024; % mol fraction; water composition
        Mw = 18.01; % g/mol; water molecular weight
        Ma = 28.97; % g/mol; air molecular weight
        Y = (yw*(Mw/Ma))/2; % kg water/kg DA; air humidity
        PartialP = Ppa*Y/((Mw/Ma)+Y); % Pa; water vapor partial pressure
        rhog = (Ppa-(0.378*PartialP))/(287.1*Tav); % kg/m3; humid air dens.

        cginlet = (yw*Patm)*1000/(R*T0); % mol/m^3; inlet humidity conc.
        tao = 0.000000000000000001; % s; time constant
        
        Cpa = 1000*(28.088+(0.197*10^(-2)*Tav)+(0.48*10^(-5)*(Tav^2))...
                -(1.965*10^(-9)*(Tav^3)))/Ma; % J/kg.K; DA heat capacity
        Cpv = 1000*(32.218+(0.192*10^(-2)*Tav)+(1.055*10^(-5)*(Tav^2))...
                -(3.593*10^(-9)*(Tav^3)))/Mw; % J/kg.K; water heat capacity
        Cpg = Cpa+(Cpv*Y); % J/kg.K; humid air heat capacity; 

        ec = 0.39; % void fraction in the column 
        ep = 0.395; % void fraction in the pellets
        rhop = 1020; % kg/m^3; solid or particle density 
        rhos = 2000; % kg/m3; pellet bulk density 
        rhow = 8238; % kg/m3; density of column wall 
        mug = (1.827*10^(-5))*((291.15+120)/(u(2)+120))*...
                ((u(2)/291.15)^(3/2)); % kg/m.s; air dvisc
        nuair = 0.00001589; % m2/s; air kinetic visc at 300K; 
        Dp = 0.0023; % m; particle diameter 
        Dc = 0.03391; % m; column internal diameter 
        Dco = 0.0381; % m; column external diameter 
        k = 0.4; % W/m.K; thermal conductivity of air at ~300 K 
        kair = 0.0263;% W/m.K; thermal conductivity of air at 300 K???
        Cps = 836.8; % J/kg.K; heat capacity of solid 
        Cpw = 473; % J/kg.K; heat capacity of the wall 
        g = 9.81; % m/s^2; gravity 
        rp = Dp/2; % m; radius of pellet 
        rc = Dc/2; % m; radius of column 
        ro = Dco/2; % m; outer radius of column with a thickness of 9mm
        dHads = -52000; % J/mol; heat of adsorption 
        Vdot = 24; % LPM; volumetric flow rate
        vg0 = Vdot*4/(PI*(Dc^2)*60000); % m/s; superficial velocity 
        vg = vg0; % m/s; assumed constant throughout from inlet
        vi = vg/ec;

        DM = 0.27; % cm2/s; molecular diffusivity 
        Dm = DM/(100^2); % m2/s; " " "
        Rep = rhog*vg*Dp/mug; % Reynolds particle number
        Sc = mug/(rhog*Dm); % Schimdt gas number
        Dz = (Dm/ec)*(20.0+0.5*Rep*Sc); % m^2/s; axial dispersion coeff
        kg = (Dm/Dp)*(2.0+1.1*(Rep^0.6)*(Sc^0.33)); % m/s; mass t.c
        Nud = 4.36; % Nusselt number in column
        hfd = k*Nud/Dc; % column overall heat t.c 
        hfp = k*Nud/Dp; % particle overall heat t.c
        Prair = 0.707; % Prandlt air number at 300K 
        Rad = (1/298.15)*(298.15-T0)*g ... 
                *(Dco^3)*Prair/(nuair^2); % Rayleight air number
        Nu = (0.60+0.387*(Rad^(1/6))*(1 ... 
                +(0.559/Prair)^(9/16))^(-8/27))^2; % Nusselt number outside
        ho = kair*Nu/(Dco); % Overall heat t.c outside
        eta = rhog*Cpg*(ec+(1-ec)*ep)+...
                (1-ec)*(1-ep)*rhos*Cps; % weighted average heat capacity
        
        % Temperature Fitted Paramters for Temperature-Dependent Toth
        a0 = 3.634*10^(-6); % mol/kg.kPa
        b0 = 2.408*10^(-7); % 1/kPa
        A = 6.852*10^3; % K
        t0 = 0.3974; % Tooth isotherm exponent coefficient
        cons = -4.199; % K
        b = b0*exp(A/u(2)); % 1/kPa
        a = a0*exp(A/u(2)); % mol/kg.kPa
        n = t0 + cons/u(2);
        cpe = u(4)*101.3/(((a*R1*u(2))^n) ... 
                -((u(4)*b*R1*u(2))^n))^(1/n); % mol/m^3; rearranged Toth

        % LDF Model Parameters
        Dso = 250; % m2/s
        Ea = -45500; % J/mol

        kads = (15*Dso*exp(Ea/(8.314*(u(2)))))/(rp^2);
        
        c = [ec; ...
            eta; ...
            ((ro^2)-(rc^2))*rhow*Cpw; ...
            1];
        % eg. f(1)=ec*Dz*DuDx(1), f(2)=k*DuDx(2), f(3)=derr*DuDx(3)~0,...
        f = [ec*Dz; ...
            k; ...
            derr; ...
            derr] .* DuDx;
        
        s = [-vg*DuDx(1)-(kads*ec)*(u(1)-cpe); ...
             -vg*rhog*Cpg*DuDx(2)-(2*hfd/rc)*(u(2)-u(3))-dHads*(kads*ec)*(u(1)-cpe); ...
             2*rc*hfd*(u(2)-u(3))-2*ro*ho*(u(3)-T0); ...
             (kads*ec/rhop)*(u(1)-cpe); ...
            ];
    end

    
    %% *1B. pdeic function definition*
    function u0 = pdeic(x)
        % MATLAB expected array output within the IC function
        u0 = [0;T0;T0;0];
    end

    
    %% *1C. pdebc function definition*
    function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,t)
        % MATLAB expected function that outputs array pl, ql, pr, and qr
        % representing the BCs in the form: 
        % p(x,t,u)+q(x,t).*f(x,t,u,Du/Dx)=0 where the l and r represent the
        % left (x=0) and right (x=H) boundaries, respectively.
        
        % left BCs with form: pl(0,t,u) + ql(0,t).*f(0,t,u,Du/Dx)=0
        pl = [vi/Dz*(cginlet-ul(1)); ... % ensures at t=0, c(t,x)=0 for all x (sharp discontinuity in graph if not)
              ul(2)-T0; ...
              ul(3)-T0; ...
              ul(4)];
        ql = [1;0;0;0];
        
        % right BCs with form: pr(H,t,u)+qr(H,t).*f(H,t,u,Du/Dx)=0
        pr = [0;0;0;0];
        qr = [1;1;1;1];
    end
    
    %% *1D. createGUI function definition*
    function createGUI
        % This function takes the results from pdepe function and contains
        % all the code relevant to creating the plots, the interactive
        % sliders, and all related aspects to the UI.
        
        fig = figure; % creates new (empty figure) window 
        sVal = 1;   % slider index, used to initialize value of indices 
                    % which represent the index value of the dependant 
                    % variable array, eg. c_g(i,j)=f(t,x), (i or j)=sindex 
                    % depending on the slider wishing to display.
        tMinutes = t/60;
        x_mm = x*1000;
        c_gRatio_xSlide = @(sVal)c_g(:,round(sVal))/cginlet;
        c_gRatio_tSlide = @(sVal)c_g(round(sVal),:)/cginlet;
        T_g_xSlide =  @(sVal)T_g(:,round(sVal));
        T_g_tSlide =  @(sVal)T_g(round(sVal),:);
        c_gExp24time = c_gExp24.(1)/60;
        c_gExp24conc = c_gExp24.(2);
        T_gExp24time = T_gExp24.(1)/60;
        T_gExp24temp = T_gExp24.(2)-273.15;
        
        %% _1D.i) createGUI section for the c/co vs t plot with x slider_
        subplot(4,2,1); % subplot function splits the figure into 4 rows, 
                        % 2 columns, and 1 is the index being set active.
                        % A 4x2 grid input represents 8 subplots within the
                        % figure, the 1st subplot is the (i,j)=(1,1)
                        % position within the grid, 2nd subplot is at the
                        % (i,j)=(1,2) position, i.e. position index is read
                        % from left to right and from up to down:
                        % e.g.  (1, 2)
                        %       (3, 4)
                        %       (5, 6)
                        %       (7, 8)
        hold on
        ploth1 = plot(tMinutes,c_gRatio_xSlide(sVal),'LineWidth',2);
        plothExp1 = plot(c_gExp24time,c_gExp24conc,'or','MarkerSize',4);
        hold off
        xlabel('Time t (min)');
        ylabel('c/c_0');
        ylim([0,1]);
        xlim([1,t_f/60]); % specifies x-axis limits in minutes
        legend('model','experimental','location','southeast');
        ax = gca; % setting gca (MATLAB keyword for current active
                  % axis) to ax for changing layout
        ax.Position(2) = ax.Position(2)-0.10;
        ax.Position(4) = ax.Position(4)+0.10;
        box on

        slider1 = uicontrol('Parent',fig,'Style','slider',...
            'Units','Normalized',...
            'Position',[ax.Position(1), ax.Position(2)-0.10,...
                        ax.Position(3), 0.03],...
            'value',sVal,...
            'min',1,'max',length(x));

        txt1 = uicontrol('Style','text',...
            'Units','Normalized',...
            'Position',[slider1.Position(1), slider1.Position(2)-0.05,...
                        slider1.Position(3), 0.03],...
            'String',sprintf('Position x = %0.1f mm',x_mm(round(sVal))));

        hLstn1 = addlistener(slider1,'ContinuousValueChange',@updateplot1);
        
        function updateplot1(~,~)
            sVal = get(slider1,'value');
            set(ploth1,'YData',c_gRatio_xSlide(sVal));
            set(txt1,'String',sprintf('Position x = %0.1f mm',...
                x_mm(round(sVal))));
        end
        
        %% _1D.ii) createGUI section for the c/co vs x plot with t slider_
        subplot(4,2,2);
        ploth2 = plot(x_mm,c_gRatio_tSlide(sVal),'LineWidth',2);
        xlabel('Position x (mm)');
        ylabel('c/c_0');
        ylim([0,1]);
        xlim([0,H*1000]);
        ax = gca; 
        ax.Position(2) = ax.Position(2)-0.10;
        ax.Position(4) = ax.Position(4)+0.10;
        box on
        
        slider2 = uicontrol('Parent',fig,'Style','slider',...
            'Units','Normalized',...
            'Position',[ax.Position(1),ax.Position(2)-0.10,...
                        ax.Position(3), 0.03],...
            'value',sVal,...
            'min',1,'max',round(length(t/60)*0.8));

        txt2 = uicontrol('Style','text',...
            'Units','Normalized',...
            'Position',[slider2.Position(1), slider2.Position(2)-0.05,...
                        slider2.Position(3), 0.03],...
            'String',sprintf('Time t = %0.1f min',tMinutes(round(sVal))));

        hLstn2 = addlistener(slider2,'ContinuousValueChange',@updateplot2);    
        
        function updateplot2(~,~)
            sVal = get(slider2,'value');
            set(ploth2,'YData',c_gRatio_tSlide(sVal));
            set(txt2,'String',sprintf('Time t = %0.1f min',...
                tMinutes(round(sVal))));
        end
        
        %% _1D.iii) createGUI section for Tg vs t plot with x slider_
        subplot(4,2,5);
        hold on
        ploth3 = plot(tMinutes,T_g_xSlide(sVal),'LineWidth',2);
        plothExp3 = plot(T_gExp24time,T_gExp24temp,'or','MarkerSize',4);
        hold off
        xlabel('Time t (min)');
        ylabel('T_g');
        yupper = round(max(max(T_g)*1.10),-1); % round mult of 10
        ylim([0,yupper]);
        xlim([1,t_f/60]);
        legend('model','experimental','location','northeast');
        ax = gca; 
        ax.Position(2) = ax.Position(2)-0.10;
        ax.Position(4) = ax.Position(4)+0.10;
        box on
        
        slider3 = uicontrol('Parent',fig,'Style','slider',...
            'Units','Normalized',...
            'Position',[ax.Position(1), ax.Position(2)-0.10,...
                        ax.Position(3), 0.03],...
            'value',sVal,...
            'min',1,'max',length(x));

        txt3 = uicontrol('Style','text',...
            'Units','Normalized',...
            'Position',[slider3.Position(1), slider3.Position(2)-0.05,...
                        slider3.Position(3), 0.03],...
            'String',sprintf('Position x = %0.1f mm',x_mm(round(sVal))));

        hLstn3 = addlistener(slider3,'ContinuousValueChange',@updateplot3);
        
        function updateplot3(~,~)
            sVal = get(slider3,'value');
            set(ploth3,'YData',T_g_xSlide(sVal));
            set(txt3,'String',sprintf('Position x = %0.1f mm',...
                x_mm(round(sVal))));
        end
        
        %% _1D.iv) createGUI section for Tg vs x plot with t slider_
        subplot(4,2,6);

        ploth4 = plot(x_mm,T_g_tSlide(sVal),'LineWidth',2);
        xlabel('Position x (mm)');
        ylabel('T_g');
        yupper = round(max(max(T_g)*1.10),-1); %round mult of 10
        ylim([0,yupper]);
        xlim([0,H*1000]);
        ax = gca; 
        ax.Position(2) = ax.Position(2)-0.10;
        ax.Position(4) = ax.Position(4)+0.10;
        box on
        
        slider4 = uicontrol('Parent',fig,'Style','slider',...
            'Units','Normalized',...
            'Position',[ax.Position(1), ax.Position(2)-0.10,...
                        ax.Position(3), 0.03],...
            'value',sVal,...
            'min',1,'max',round(length(tMinutes)*0.8));

        txt4 = uicontrol('Style','text',...
            'Units','Normalized',...
            'Position',[slider4.Position(1), slider4.Position(2)-0.05,...
                        slider4.Position(3), 0.03],...
            'String',sprintf('Time t = %0.1f min',t(round(sVal))/60));

        hLstn4 = addlistener(slider4,'ContinuousValueChange',@updateplot4);    
        
        function updateplot4(~,~)
            sVal = get(slider4,'value');
            set(ploth4,'YData',T_g_tSlide(sVal));
            set(txt4,'String',sprintf('Time t = %0.1f min',...
                tMinutes(round(sVal))));
        end
    end
end
