% ME140 PROJECT 4: FUEL CELL
% --------------------------
% Jon Renslo, Emily Bohl, Frankie Willcox, Natasha Berk, Kendall Fagan

% ----------------------------------------------
% Part 2b: Efficiency as a Function of Pressure
% ----------------------------------------------
close all; clear; clc;

% CONSTANTS
Pmax = 40;                   % [atm]                      

% Unit Conversions
C_TO_K = 273.15;
ATM_TO_PA = 101325; 

% Molecular weights             [g/mol]
MW_O2 = 2*15.9994;
MW_N2 = 2*14.0067; 
MW_H2 = 2*1.00794; 
MW_H2O = MW_H2 + MW_O2*.5; 

% Ethalpies of formation        [J/mol]
hf_O2 = 0;  
hf_H2 = 0; 
hf_N2 = 0; 
hf_H2OVAP = -241820; 
hf_H2OLIQ = -285830;  

% Entropies of Formation        [J/mol*K]
sf_H2OVAP = 188.83; 
sf_H2OLIQ = 69.92; 
sf_O2 = 205.04;
sf_N2 = 191.61;
sf_H2 = 130.68;

% Heating Values                [J/g]
LHV = 120*10^3; 
HHV = 141.8*10^3; 

% Atmospheric Conditions
Tref = 298;                  % [K]
Patm = 1 * ATM_TO_PA;


T_range = [80 220 650 800] + C_TO_K;

% Constant Cp of Liquid Water
Cp_H2OLIQ = 75.96; 
cp_const_T = @(T) (Cp_H2OLIQ./T); 

for j = 1:4
    T = T_range(j);
    for i = 1:Pmax                   
        lambda = 2;                 % ASSUME: 100% saturated
        
        % Increment Pressure
        Pmix = i*ATM_TO_PA;           % [Pa]
        
        % Initialize alpha & beta
        alpha = 0; 
        beta = 1 + alpha; 
        
        % Mols of Each Species (found from balanced chemical equation) 
        % ------------------------------------------------------------
        % Same for Reactants & Products (N2)
        mol_n2_react = 0.5*lambda*3.76; % stays same for reactants & products 
        mol_n2_prod = mol_n2_react;
        
        % Reactants (H2,O2,H2O_liquid) 
        mol_h2_react = 1; % h2 = fuel
        mol_o2_react = 0.5*lambda; 
        mol_h2oliq_react = alpha; 
        
        % Products ( O2,N2(defined above),H2O_vapor(defined as beta below) )
        mol_o2_prod = 0.5*(lambda-1);  
  
        % Mass of Fuel (H2)
        % -----------------
        m_fuel = mol_h2_react*MW_H2;
        
        beta(i) = 1+alpha; 
        Psat(i) = PsatW(T); 
        Pv(i) = Pmix*(beta(i))/(beta(i)+mol_o2_prod + mol_n2_prod); 
        if(Pv(i) <= Psat(i)) 
            beta(i) = alpha+1;
            gamma(i) = alpha-beta(i)+1;
        else 
            beta(i) = (mol_o2_prod + mol_n2_prod)./(Pmix./Psat(i)-1); 
            gamma(i) = alpha -beta(i) +1;
        end 
        
        prod_Ntotal = beta(i) + gamma(i) + mol_n2_prod + mol_o2_prod; 
       
       %----------------------------------------------
       % STEP 1: Enthalpy Change (dh) for each Species
       %----------------------------------------------
       H2_dH(i) = hf_H2 + dh(Tref,T,'H2');
       O2_dH(i) = hf_O2 + dh(Tref,T,'O2'); 
       N2_dH(i) = hf_N2 + dh(Tref,T,'N2');
       H2O_gas_dH(i) = hf_H2OVAP + dh(Tref,T,'H2Ovap');
       H2O_liquid_dH(i) = hf_H2OLIQ + Cp_H2OLIQ*(T-Tref);
      
       %-----------------------------------------------------------------
       % STEP 2: Partial Pressure (Ppartial), Entropy Change (dS), Gibbs
       % Free Energy (dG)
       %-----------------------------------------------------------------
       
       % REACTANTS (Note: H2 is separated by membrane, does not contribute
       % to total pressure--> don't include in denominator when calculating partial pressures)
       % ---------
       % Partial Pressure
       Ppartial_O2_reactant(i) = Pmix*(mol_o2_react/(mol_o2_react+mol_n2_react));
       Ppartial_N2_reactant(i) = Pmix*(mol_n2_react/(mol_o2_react+mol_n2_react));
       
       % dS
       dS_h2_react(i) = sf_H2           + dS(Tref,T,Pmix                   ,Patm,'H2'); 
       dS_n2_react(i) = sf_N2           + dS(Tref,T,Ppartial_N2_reactant(i),Patm,'N2');
       dS_o2_react(i) = sf_O2           + dS(Tref,T,Ppartial_O2_reactant(i),Patm,'O2');
       dS_h2o_react(i) = sf_H2OVAP      + dS(Tref,T,Pmix                   ,Patm,'H2Ovap');
       % This is different line here -JR
       % dG
       dG_h2_react(i) =     H2_dH(i)      - T*dS_h2_react(i);
       dG_o2_react(i) =     O2_dH(i)      - T*dS_o2_react(i);
       dG_n2_react(i) =     N2_dH(i)      - T*dS_n2_react(i);
       dG_h2ovap_react(i) = H2O_gas_dH(i) - T*dS_h2o_react(i);
       
       % PRODUCTS
       % ---------
       % Partial Pressure
       Ppartial_O2_product(i) = Pmix*(mol_o2_prod /(mol_o2_prod + mol_n2_prod+beta(i)));
       Ppartial_N2_product(i) = Pmix*(mol_n2_prod /(mol_o2_prod + mol_n2_prod+beta(i)));
       Ppartial_H2O_product(i) = Pmix*(beta(i)    /(mol_o2_prod + mol_n2_prod+beta(i)));
       
       % dS
       dS_o2_prod(i) = sf_O2           + dS(Tref,T,Ppartial_O2_product(i),Patm,'O2');
       dS_n2_prod(i) = sf_N2           + dS(Tref,T,Ppartial_N2_product(i),Patm,'N2');
       dS_h2ovap_prod(i) = sf_H2OVAP + dS(Tref,T,Ppartial_H2O_product(i),Patm,'H2Ovap');
       dS_h2oliq_prod(i) = integral(cp_const_T,Tref,T) + sf_H2OLIQ;
       
       % dG
       dG_o2_prod(i) = O2_dH(i) - T*dS_o2_prod(i);
       dG_n2_prod(i) = N2_dH(i) - T*dS_n2_prod(i);
       dG_h2ovap_prod(i) = H2O_gas_dH(i) - T*dS_h2ovap_prod(i);
       dG_h2oliq_prod(i) = H2O_liquid_dH(i) - T*dS_h2oliq_prod(i);   
       
       G_react(i) = mol_h2_react*dG_h2_react(i) + mol_o2_react*dG_o2_react(i) + mol_n2_react*dG_n2_react(i) + alpha*dG_h2ovap_react(i);
       G_prod(i) = mol_o2_prod*dG_o2_prod(i) + mol_n2_prod*dG_n2_prod(i) + beta(i)*dG_h2ovap_prod(i) + gamma(i)*dG_h2oliq_prod(i);
       dG_rxn(i) = G_prod(i) - G_react(i);
       
       eta_HHV(i) = -(dG_rxn(i))./(m_fuel*HHV);
       eta_LHV(i) = -(dG_rxn(i))./(m_fuel*LHV);
       eta_carnot(i) = (T-Tref)./T;
    end
    eta_press(j,:) = eta_LHV;
end 

P_range = linspace(1,Pmax,Pmax) * ATM_TO_PA; % [Pa]
f = figure
hold on
plot(P_range,eta_press(1,:),'r',P_range,eta_press(2,:),'y',P_range,eta_press(3,:),'g',P_range,eta_press(4,:),'b')
legend('T = 80šC','T = 220 C','T = 650 šC','T = 800 šC','location','southeast')
title('Part 2b: Efficiency on a LHV basis as a Function of Pressure')
xlabel('Pressure [Pa]')
ylabel('Efficiency \eta_{LHV}')
plotfixer
grid on
saveas(f,'../plots/Part2-2','jpeg');