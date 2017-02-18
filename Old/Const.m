function y = Const(str)

    switch str
        case 'R'
            y = 8.3144621;      % J/mol.K
        case 'k'
            y = 1.3806488e-23;  % m2.kg/s2.K
        case 'pi'
            y = pi();
        case 'g'
            y = 9.81;           % m/s^2; gravitational constant
    end
end