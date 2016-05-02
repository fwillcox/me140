% ME 140 Project #5
% FUEL CELL EVALUATION & HYRDOGEN PRODUCTION ANALYSIS
% Frankie Willcox, Jon Renslo, Kendall Fagan, Emily Bohl, Natasha Berke

% ASSUME:
% (i)  mol_H2 = 1

% Constants
G_TO_KG = 10^3;

% Molar Masses
MM_h = 1.00794; % g/mol
MM_o = 15.9994;
MM_n = 14.0067;
MM_h2o = 2*MM_h + MM_o;
MM_air = 28.97;

global PERMIN_TO_PERSEC PERHR_TO_PERSEC G_PER_KG LHV F N_TO_O SCF_TO_MOLS C_TO_K
defineGlobals();
mol_H2 = 1;
savePlots = 1;

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
T1 = T1_C + C_TO_K;                                                % [K],  T1, air into stack
T2_C =     [42.5 45.8 45.8 48.4 50.3 51.9 53.3 53.9 54.3];         
T2 = T2_C + C_TO_K;                                                % [K],  T2, air out of stack 
Tstack_C = [40.7 41.3 42.5 42.9 44.6 45.6 46.9 46.9 47.6];  
Tstack = Tstack_C + C_TO_K;

% NOTE: T3-T5 are not needed for now
% T3_C =     [48.0 47.1 48.6 48.9 50.4 51.1 51.2 51.1 51.1];       % T3, water reservoir DON'T USE!
% T3 = T3_C + C_TO_K;
% T4_C =     [48.0 47.2 48.2 48.9 50.4 51.1 51.2 51.1 51.1];       % T4, water into stack
% T4 = T4_C + C_TO_K;
% T5_C =     [40.7 41.3 42.5 42.9 44.6 45.6 46.9 46.9 47.6];       % T5, water into heat exchanger
% T5 = T5_C + C_TO_K;

% Mass Flow Rates (TODO: check what units the mdots should be in)
mdot_h2_scf =  [2.50 6.20 10.5 14.3 18.2 22.0 24.6 25.0 26.1];     % scf/hr (standard cubic feet/hour)
mdot_air_scf = [0.75 1.10 1.45 1.81 2.55 3.10 3.30 3.25 3.40];     % scf/min

mdot_h2 = mdot_h2_scf   * SCF_TO_MOLS * PERHR_TO_PERSEC * MM_h2  / G_TO_KG;  % kg/s	
mdot_air = mdot_air_scf * SCF_TO_MOLS * PERMIN_TO_PERSEC* MM_air / G_TO_KG;  % kg/s

mdot_h2o = 40 /G_TO_KG;                                            % kg/s (ASSUMPTION for Part 3)


% CALCULATIONS
% ------------
% Power (USE: p = i*v) 
% NOTE: "Accessories" include H2O pump, air pump, H2 vent, & controller                 
p_load =  i_load  .* v_load;                                        % a.k.a. "load" (power delivered to resistor bank)
p_stack = i_stack .* v_stack;
p_access = p_stack - p_load;                                        % Acessory Power, i.e. power used to run controls. Pstack-Pload

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
plot(p_load,mdot_h2,p_load,mdot_air);
title('Mass Flow Rate as a Function of Load');
xlabel('Load []'); ylabel('Mass Flow Rate []');
legend('mdot_{H}','mdot_{air}'); plotfixer();grid on

% ---------------------------
% Part 2: Reduced-Data Plots
% ---------------------------
% SOURCE: LEC 8, SLIDES 21 & 22

% 1st & 2nd Law Efficiencies (eta_I & eta_II) & Inefficiencies (Idot)
% Stack
[etaI_stack ,etaII_stack, Idot_stack] = findEtas(Tstack, p_stack);

% Entire System (Load)
[etaI_load ,etaII_load, Idot_load] =    findEtas(Tstack, p_load);

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




