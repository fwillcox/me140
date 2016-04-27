function [ dS ] = dS( T1,T2,P1,P2,gas )
% finds change in entropy between 2 temperatures for a variable Cp gas

R = 8.31441;    % J/(K*mol)
p = [0 0 0 0];

   % Coefficients of Cp=Cp(T) for each gas
   switch(gas) 
    case 'air'
        p = [-1.966*10^-9 0.4802*10^-5 0.1967*10^-2 28.11]; 

    case 'H2Ovap' 
        p = [-3.595*10^-9 1.055*10^-5 0.1923*10^-2 32.24]; 

    case 'N2'
        p = [-2.873*10^-9 0.8081*10^-5 -0.1571*10^-2 28.90]; 

    case 'O2'
        p = [1.312*10^-9 -0.7155*10^-5 1.520*10^-2 25.48]; 
    
    case 'H2'
        p = [-0.8704*10^-9 0.4003*10^-5 -0.1916*10^-2 29.11]; 
        
   end
   
   % Entropy Equation: dS = integral from T1 to T2 of(cp(T)/T) - R*ln(P1/P2)
    cp_var_by_T = @(T) ((p(1).*T.^2 + p(2).*T + p(3) + p(4)./T)); %[kJ/kmol*K]
    dS = integral(cp_var_by_T,T1,T2) - R*log(P1/P2);
end

