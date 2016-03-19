function dwdt = LDF(type)
Ru = RGAS('J');
ISO = Isotherm(type);

if strcmp(type,'ZEO13X')
    w_eq = ISO.w_eq;
    k1   = 4.004e-2; % 1/s
    k2   = 905.8; % K

    dwdt = k1.*exp(-k2./T).*(w_eq - w);

elseif strcmp(type,'SILGEL')
    w_eq = ISO.w_eq;
    Dso  = 2.54e-4; % m2/s
    Rp   = 1.74e-4; % m
    Ea   = 42000; % J/mol

    dwdt = 15.*Dso./Rp^2.*exp(-Ea./(Ru.*T)).*(w_eq - w);
end
