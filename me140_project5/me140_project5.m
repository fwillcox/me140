% ME 140 Project #5
% FUEL CELL EVALUATION & HYRDOGEN PRODUCTION ANALYSIS
% Frankie Willcox, Jon Renslo, Kendall Fagan, Emily Bohl, Natasha Berke

% ASSUME:
% (i)  

% Constants
PERMIN_TO_PERHR = 60;
LHV = 120.0*10^6;                                                   % J/kg,  Lower Heating Value 
F = 96485;                                                          % C/(mol of e-), Faraday's Constant
mol_H2 = 1;

% --------------------------------
% Part 1: Raw Data Plots vs. Load
% --------------------------------
% Currents (load & stack)
i_load =  [0.00 15.06 27.25 36.48 45.1 52.1 56.3 57.6 56.4];        % Amps
i_stack = [4.82 21.40 35.65 47.20 59.8 69.7 77.0 79.0 80.0];

% Potentials (load & stack)
v_load =  [17.07 15.05 14.08 13.10 12.07 11.27 10.31 9.87  9.05 ];  % Volts
v_stack = [17.09 15.22 14.26 12.98 12.42 11.60 10.73 10.21 9.48];

% Temperatures from Thermocouple Readings [C]
% KEY: (Kendall please fill in with photo you took)
T1 =     [42.8 42.9 46.1 48.5 50.5 52.8 54.8 55.8 56.5];            % T1, air into stack
T2 =     [42.5 45.8 45.8 48.4 50.3 51.9 53.3 53.9 54.3];            % T2, air out of stack
Tstack = [40.7 41.3 42.5 42.9 44.6 45.6 46.9 46.9 47.6];            % **USE Tstack for Psat

% NOTE: T3-T5 are not needed for now
%T3 =     [48.0 47.1 48.6 48.9 50.4 51.1 51.2 51.1 51.1];           % T3, water reservoir DON'T USE!
%T4 =     [48.0 47.2 48.2 48.9 50.4 51.1 51.2 51.1 51.1];            % T4, water into stack
%T5 =     [40.7 41.3 42.5 42.9 44.6 45.6 46.9 46.9 47.6];            % T5, water into heat exchanger


% CALCULATIONS
% ------------
% Power (USE: p = i*v) 
% NOTE: "Accessories" include H2O pump, air pump, H2 vent, & controller                 
p_load =  i_load  .* v_load;                                        % a.k.a. "load" (power delivered to resistor bank)
p_stack = i_stack .* v_stack;
p_access = p_stack - p_load;                                        % Acessory Power, i.e. power used to run controls. Pstack-Pload

figure(1)
plot(p_load,i_load,'c',p_load,i_stack,'bp--');
title('Current as a Function of Load');
xlabel('Load [Watts]'); ylabel('Current [Amps]');
legend('I_{load}','I_{stack}'); plotfixer(); grid on

figure(2) 
plot(p_load,v_load,'c',p_load,v_stack,'bp--');
title('Potential as a Function of Load');
xlabel('Load [Watts]'); ylabel('Potential [V]');
legend('V_{load}','V_{stack}'); plotfixer();grid on

figure(3)
plot(p_load,p_stack,'c',p_load,p_access,'bp--');
title('Stack and Accessory Power as a Function of Load');
xlabel('Load [Watts]'); ylabel('Power [Watts]');
legend('P_{stack}','P_{accessory}'); plotfixer();grid on

figure(4)
plot(p_load,mdot_h2,'c',p_load,mdot_air,'bp--');
title('Mass Flow Rate as a Function of Load');
xlabel('Load [Watts]'); ylabel('Mass Flow Rate []');
legend('mdot_{H}','mdot_{air}'); plotfixer();grid on

% ---------------------------
% Part 2: Reduced-Data Plots
% ---------------------------
% SOURCE: LEC 8, SLIDES 21 & 22

% 1st & 2nd Law Efficiencies (eta_I & eta_II) & Inefficiencies (Idot)
% Stack
[etaI_stack ,etaII_stack, Idot_stack] = findEtas(Tstack);

% Entire System (Load)
[etaI_load ,etaII_load, Idot_load] = findEtas(Tstack);

figure(5)
plot(p_load,etaI_stack,'c',p_load,etaI_load,'bp--');
title('First Law Efficiency as a Function of Load');
xlabel('Load [Watts]'); ylabel('Efficiency, eta_{I}');
legend('eta_{I,stack}','eta_{I,system}'); plotfixer(); grid on

figure(6) 
plot(p_load,etaII_stack,'c',p_load,etaII_load,'bp--');
title('Second Law Efficiency as a Function of Load');
xlabel('Load [Watts]'); ylabel('Efficiency, eta_{II}');
legend('eta_{II,stack}','eta_{II,system}'); plotfixer(); grid on

figure(7)
plot(p_load,p_stack,'c',p_load,p_load,'bp--');
title('Power Loss/Inefficiences as a Function of Load');
xlabel('Load [Watts]'); ylabel('Power Loss/Inefficiencies, Idot [Watts]');
legend('Idot_{stack}','Idot_{system}'); plotfixer();grid on



