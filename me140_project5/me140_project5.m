% ME 140 Project #5
% FUEL CELL EVALUATION & HYRDOGEN PRODUCTION ANALYSIS
% Frankie Willcox, Jon Renslo, Kendall Fagan, Emily Bohl, Natasha Berke

% ASSUME:
% (i)  mol_H2 = 1

global PERMIN_TO_PERSEC PERHR_TO_PERSEC G_PER_KG LHV F N_TO_O SCF_TO_MOLS ...
    C_TO_K PSI_TO_PA MM_h MM_h2 MM_o MM_n MM_h2o MM_air PATM
defineGlobals();
mol_H2 = 1;
savePlots = 0;

% --------------------------------
% Part 1: Raw Data Plots vs. Load
% --------------------------------
% Currents (load & stack)
i_load =  [0.00 15.06 27.25 36.48 45.1 52.1 56.3 57.6 56.4];       % [Amps]
i_stack = [4.82 21.40 35.65 47.20 59.8 69.7 77.0 79.0 80.0];


% Potentials (load & stack)
v_load =  [17.07 15.05 14.08 13.10 12.07 11.27 10.31 9.87  9.05 ]; % [Volts]
v_stack = [17.09 15.22 14.26 12.98 12.42 11.60 10.73 10.21 9.48];


% Temperatures from Thermocouple Readings [C]
% KEY: (Kendall please fill in with photo you took)
T1_C = [42.8 42.9 46.1 48.5 50.5 52.8 54.8 55.8 56.5];
T2_C =     [42.5 45.8 45.8 48.4 50.3 51.9 53.3 53.9 54.3];
Tstack_C = [40.7 41.3 42.5 42.9 44.6 45.6 46.9 46.9 47.6];  

T1 = T1_C + C_TO_K;                                                % [K],  T1, air into stack
T2 = T2_C + C_TO_K;                                                % [K],  T2, air out of stack 
Tstack = Tstack_C + C_TO_K;                                        % [K],  metal plates on the stack

% NOTE: T3-T5 are not needed for now
% T3_C =     [48.0 47.1 48.6 48.9 50.4 51.1 51.2 51.1 51.1];       % T3, water reservoir DON'T USE!
% T3 = T3_C + C_TO_K;
% T4_C =     [48.0 47.2 48.2 48.9 50.4 51.1 51.2 51.1 51.1];       % T4, water into stack
% T4 = T4_C + C_TO_K;
% T5_C =     [40.7 41.3 42.5 42.9 44.6 45.6 46.9 46.9 47.6];       % T5, water into heat exchanger
% T5 = T5_C + C_TO_K;


% Mass Flow Rates (TODO: check what units the mdots should be in)
mdot_total_scf = [0.75 1.10 1.45 1.81 2.55 3.10 3.30 3.25 3.40];   % scf/min
mdot_fuel_scf =  [2.50 6.20 10.5 14.3 18.2 22.0 24.6 25.0 26.1];   % scf/hr (standard cubic feet/hour)

mdot_total = mdot_total_scf * SCF_TO_MOLS * PERMIN_TO_PERSEC * ( MM_air / G_PER_KG);  % kg/s
mdot_fuel = mdot_fuel_scf * SCF_TO_MOLS * PERHR_TO_PERSEC * (MM_h2  / G_PER_KG);  % kg/s	
mdot_h2o = 40 /G_PER_KG;                                                        % kg/s


% Pressures
Pfuel_psi =  [2.9 2.9 3.1 3.3 3.30 3.20 3.00 3.0 3.1];              % [psi] (gauge)
Ptotal_psi = [0.2 0.3 0.6 0.7 1.15 1.25 1.35 1.3 1.5];              % [psi] (gauge), pressure of combined air and H2O after humidifier
Pfuel =  Pfuel_psi  .* PSI_TO_PA + PATM;                            % [Pa] 
Ptotal = Ptotal_psi .* PSI_TO_PA + PATM;                            % [Pa] 

% CALCULATIONS
% ------------
% Power (USE: p = i*v) 
% NOTE: "Accessories" include H2O pump, air pump, H2 vent, & controller                 
p_load =  i_load  .* v_load;                                        % [W] = [kg*m^2*s^-3], a.k.a. "load" (power delivered to resistor bank)
p_stack = i_stack .* v_stack;
p_access = p_stack - p_load;                                        % [W], Acessory Power, i.e. power used to run controls. Pstack-Pload

f1 = figure(1);
plot(p_load,i_load,p_load,i_stack);
title('Current as a Function of Load');
xlabel('Load []'); ylabel('Current []');
legend('I_{load}','I_{stack}'); plotfixer(); grid on;

f2 = figure(2); 
plot(p_load,v_load,p_load,v_stack);
title('Potential as a Function of Load');
xlabel('Load []'); ylabel('Potential []');
legend('V_{load}','V_{stack}'); plotfixer();grid on;

f3 = figure(3);
plot(p_load,p_stack,p_load,p_access);
title('Stack and Accessory Power as a Function of Load');
xlabel('Load []'); ylabel('Power []');
legend('P_{stack}','P_{accessory}'); plotfixer();grid on;

f4 = figure(4);
plot(p_load,mdot_fuel,p_load,mdot_total);
title('Mass Flow Rate as a Function of Load');
xlabel('Load []'); ylabel('Mass Flow Rate []');
legend('mdot_{H}','mdot_{air}'); plotfixer();grid on

% ---------------------------
% Part 2: Reduced-Data Plots
% ---------------------------
% SOURCE: LEC 8, SLIDES 21 & 22

% 1st & 2nd Law Efficiencies (eta_I & eta_II) & Inefficiencies (Idot)
% Stack
[etaI_stack ,etaII_stack, Idot_stack] = findEtas(mdot_total, mdot_fuel, Ptotal, Pfuel, T2, p_stack);

% Entire System (Load)
[etaI_load ,etaII_load, Idot_load] =    findEtas(mdot_total, mdot_fuel, Ptotal, Pfuel, T2, p_load);

f5 = figure(5);
plot(p_load,etaI_stack,'c',p_load,etaI_load,'bp--');
title('First Law Efficiency as a Function of Load');
xlabel('Load [Watts]'); ylabel('Efficiency, eta_{I}');
legend('eta_{I,stack}','eta_{I,system}'); plotfixer(); grid on;


f6 = figure(6); 
plot(p_load,etaII_stack,'c',p_load,etaII_load,'bp--');
title('Second Law Efficiency as a Function of Load');
xlabel('Load [Watts]'); ylabel('Efficiency, eta_{II}');
legend('eta_{II,stack}','eta_{II,system}'); plotfixer(); grid on;


f7 = figure(7);
plot(p_load,p_stack,'c',p_load,p_load,'bp--');
title('Power Loss/Inefficiences as a Function of Load');
xlabel('Load [Watts]'); ylabel('Power Loss/Inefficiencies, Idot [Watts]');
legend('Idot_{stack}','Idot_{system}'); plotfixer();grid on;


if(savePlots ==1) 
    saveas(f1,'../plots5/1-CurrentbyLoad','png');
    saveas(f2,'../plots5/2-VbyLoad','png'); 
    saveas(f3,'../plots5/3-PowerbyLoad','png'); 
    saveas(f4,'../plots5/4-massbyload','png'); 
    saveas(f5,'../plots5/5-FirstLaw','png'); 
    saveas(f6,'../plots5/6-SecondLaw','png'); 
    saveas(f7,'../plots5/7-PowerLoss','png'); 
end


%% Part A, Section 3 - Emily
% Comparing First Law Efficiencies

% Typical modern Diesel engine = 42% (chose diesel truck because it's
% better than a car and worse than a freight ship)
% Source: Slide 3, http://www.sae.org/events/gim/presentations/2011/RolandGravel.pdf
dieselEff = 0.42;

% Typical gasoline hybrid engine = max of 40%
% Source: Toyota Hybrid Vehicles, http://www.toyota-global.com/innovation/environmental_technology/hybrid/
hybridEff = 0.40;

% Overall First Law Efficiency of the PEM Fuel Cell = Stack Efficiency
% Plotting to compare
figure(8)
plot(p_load, etaI_stack, 'c', p_load, dieselEff, 'bp--', p_load, hybridEff, 'gd');
title('Comparing 1st Law Efficiency: PEM Fuel Cell, Diesel, and Gasoline Hybrid');
xlabel('Load [Watts]'); ylabel('Efficiency, eta_{I}');
legend('eta_{I,stack}','eta_{I,Diesel}', 'eta_{I,Hybrid}'); plotfixer(); grid on;


% TODO: FIGURE OUT SCALE-UP FOR TOYOTA HYBRID
% TODO: COMMENT ON ACCESSORY/FUEL SYSTEMS REQUIRED FOR THAT SCALE UP

%% Part B

% Part B, Section 1 - Emily & Kendall
% Calculating Kp Values
% Formulas from https://coursework.stanford.edu/access/content/group/Sp16-ME-140-01/Lecture%20Slides/Lecture%2013.pdf

% SMR: CH4 + H2O --> CO + 3H2
% v values are stoichiometric coefficients
v_CO_SMR = 1;
v_H2_SMR = 3;
v_H2O_SMR = 1;
v_CH4_SMR = 1;

% Calculating Kp for SMR
% Nv_CO = mm
% SMRnumKp = 




% WGS: H2O + CO --> H2 + CO2
v_H2_WGS = 1;
v_CO2_WGS = 1;
v_H2O_WGS = 1;
v_CO_WGS = 1;

T_B1 = linspace(25, 1200, 100); %Temperature for part B1 = T_B1
T_B1 = T_B1 + C_TO_K;
%NOTE: Stanford pressure, is usually defined as 100,000, however in energyF
%we have standard presssure as 101300. Because the pressure needs to
%cancel out, I have changed this pressure to 101300, however, we should
%perhaps consider changing the reference pressure in energyF to 100,000Pa.
P_ref = 101300; %This is the pressure defined for standard conditions. Standard conditions are what we need because that is what the little zero indicates in the equation for g. 
R_u = 8.314; %Universal gas constant

%G_reaction = G_products - G_reactants
g_SMR = (gEng(T_B1, P_ref, 'co',v_CO_SMR) + gEng(T_B1, P_ref, 'h2',v_H2_SMR)) - ...
    (gEng(T_B1, P_ref, 'h2ovap',v_H2O_SMR) + gEng(T_B1, P_ref, 'ch4',v_CH4_SMR));

g_WGS = (gEng(T_B1, P_ref, 'h2',v_H2_WGS) + gEng(T_B1, P_ref, 'co2',v_CO2_WGS)) - ...
    (gEng(T_B1, P_ref, 'h2ovap',v_H2O_WGS) + gEng(T_B1, P_ref, 'co',v_CO_WGS));


kp_SMR = exp(-g_SMR ./ (R_u .* T_B1)); %increases with temp
kp_WGS = exp(-g_WGS ./ (R_u .* T_B1)); %decrease with temp


%Prep for plot
%convert back to celcius
T_B1 = T_B1 - C_TO_K;
%find index of where kp=10^-3 and kp = 10^3, as the problem asks that we
%limit the graph to this range
[~,i_min_SMR] = min(abs(kp_SMR - 10^-3));
[~,i_max_SMR] = min(abs(kp_SMR - 10^3)); %yes, this is supposed to use min() to find the max ;P
[~,i_min_WGS] = min(abs(kp_WGS - 10^3));
[~,i_max_WGS] = min(abs(kp_WGS - 10^-3));

%Plot
figure(9)
semilogy(T_B1(i_min_SMR:i_max_SMR), kp_SMR(i_min_SMR:i_max_SMR), ...
    T_B1(i_min_WGS:i_max_WGS), kp_WGS(i_min_WGS:i_max_WGS));
xlabel('Temperature [C]')
ylabel('Equilibrium Constant')
legend('SMR', 'WGS')
title('Part B.1: Equilibrium Constant vs. Temperature')
