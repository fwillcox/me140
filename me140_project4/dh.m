function [ dh ] = dh( T1,T2,gas )
% Finds change in enthalpy of a gas from two input temperatures

p = [0 0 0 0];
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
    
cp_var = @(T) ((p(1).*T.^3 + p(2).*T.^2 + p(3).*T + p(4))); %[J/kmol*K]
dh = integral(cp_var,T1,T2); 
end

