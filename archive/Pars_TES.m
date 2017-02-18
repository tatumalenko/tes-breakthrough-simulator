function Sol = Pars_TES
    Sol.R  = Const('R'); % J/mol.K or Pa.m^3/mol.K
    Sol.PI = Const('pi');
    Sol.kB = Const('k');
    Sol.g  = Const('g');
    
    % BED CHARACTERISTICS
    Sol.H  = 0.0695;        % m; height of column 
    % Sol.eb = 0.4 + 0.05*(Sol.Dp/Sol.Dint) + 0.412*(Sol.Dp/Sol.Dint)^2; % Sol.Dp/Sol.Dint < 0.5 (Dixon 1988)
    Sol.eb    = 0.39;       % --; bed void fraction = 0.375 + 0.34*Dp/Dint (Jeshar eqn)
    Sol.Dint  = 0.03391;    % m; column internal diameter 
    Rint = Dint/2;          % m; radius of column 
    Sol.Dext  = 0.0381;     % m; column external diameter 
    Rext = Dext/2;          % m; outer radius of column with a thickness of 9mm
    Sol.rhob  = 900;        % (1-eb)*rhop where rhop = (1-ep)*rhos
    
    % INLET FLOW CHARACTERISTICS
    Sol.T0      = 298.15;                   % K; initial temperature,% K; ambient temperature
    Sol.Tginlet = Sol.T0;
    Sol.Patm    = 1.1;                      % atm; operating pressure, % Pa; atmospheric pressure (outside column)
    Sol.Ppa     = Sol.Patm*101325;          % Pa; operating Pressure 
    Sol.yw      = 0.024;                    % mol fraction; water composition
    Sol.cginlet = (Sol.yw*Sol.Patm)*1000/(0.082*Sol.T0);   % mol/m^3; inlet humidity conc.
    Sol.Mw      = 18.01;                    % g/mol; water molecular weight
    Sol.Ma      = 28.97;                    % g/mol; air molecular weight
    Sol.Vdot    = 24;                       % LPM; volumetric flow rate 
    Sol.vg      = Sol.Vdot*4/(Sol.PI*(Sol.Dint^2)*60000); % m/s; bulk interstitial velocity
    Sol.dHads   = 52000;                    % J/mol; heat of adsorption 
    Sol.dipmo_v = 1.85;
    Sol.dipmo_a = 0;
    Sol.Vlmnb_v = 5.3198e4;            % mol/m^3; 
    Sol.Vlmnb_a = 3.0215e4;            % mol/m^3;
    Sol.Tnb_v   = 373.1243;
    Sol.Tnb_a   = 78.9030;     
%     Sol.cginlet = 2.5;
%     Sol.rhog    = 1000;
%     Sol.Dz      = 5.9e-20;
%     Sol.De      = 2.4e-11;
%     Sol.kf      = 8.2e-10;
%     Sol.vDot    = 4;
%     Sol.vg      = 4e-6/60/Sol.PI*Sol.vDot/Sol.Dint^2;
%     Sol.hf      = 8.217;
%     Sol.kg      = 0.026;
%     Sol.Cpg     = 1005;
%     Sol.lambda  = Sol.kg;
%     Sol.hw      = 8.217;
    
    % PELLET CHARACTERISTICS
    Sol.Dp   = 0.0023;   % m; mean particle diameter
    Sol.Rp   = Sol.Dp/2;
    Sol.rp   = 1.75e-9;  % m; mean pore radius
    Sol.ep   = 0.53;     % --; particle void fraction
    Sol.rhos = 1970;     % kg/m3; particle density
    Sol.Cps  = 836.8;    % J/kg.K; heat capacity of solid 
    Sol.ks   = 0.147;    % W/m.K
%     Sol.hf = 30.77;
    Sol.rhop = (1 - Sol.ep)*Sol.rhos;

    % COLUMN WALL CHARACTERISTICS
    Sol.rhow = 8238;  % kg/m3; density of column wall 
    Sol.Cpw  = 473;   % J/kg.K; heat capacity of the wall 
    Sol.kw   = 17;    % W/m.K; thermal conductivity of wall
    % Sol.Uo   = 2.365;
end