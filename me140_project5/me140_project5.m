% ME 140 Project #5
% FUEL CELL EVALUATION & HYRDOGEN PRODUCTION ANALYSIS
% Frankie Willcox, Jon Renslo, Kendall Fagan, Emily Bohl, Natasha Berke

% Constants
LHV = 120.0*10^6;                    % J/kg,  Lower Heating Value  

% --------------------------------
% Part 1: Raw Data Plots vs. Load
% --------------------------------
% Currents (load & stack)
i_load =  [0.00 15.06 27.25 36.48 45.1 52.1 56.3 57.6 56.4];        % Amps
i_stack = [4.82 21.40 35.65 47.20 59.8 69.7 77.0 79.0 80.0];

% Potentials (load & stack)
v_load =  [17.07 15.05 14.08 13.10 12.07 11.27 10.31 9.87  9.05 ];  % Volts
v_stack = [17.09 15.22 14.26 12.98 12.42 11.60 10.73 10.21 9.48];

% Mass Flow Rates (H & air)
mdot_H2 =  [2.50 6.20 10.5 14.3 18.2 22.0 24.6 25.0 26.1];          % scf/h (standard cubic feet/hour)
mdot_air = [0.75 1.10 1.45 1.81 2.55 3.10 3.30 3.25 3.40];          % scf/m
mdot_cool = 40;                                                     % g/s (ASSUMPTION for Part 3)

% Temperatures from Thermocouple Readings [C]
% KEY: (Kendall please fill in with photo you took)
T1 =     [42.8 42.9 46.1 48.5 50.5 52.8 54.8 55.8 56.5]; % T1, air into stack
T2 =     [42.5 45.8 45.8 48.4 50.3 51.9 53.3 53.9 54.3]; % T2, air out of stack
T3 =     [48.0 47.1 48.6 48.9 50.4 51.1 51.2 51.1 51.1]; % T3, water reservoir
T4 =     [48.0 47.2 48.2 48.9 50.4 51.1 51.2 51.1 51.1]; % T4, water into stack
T5 =     [40.7 41.3 42.5 42.9 44.6 45.6 46.9 46.9 47.6]; % T5, water into heat exchanger
Tstack = [40.7 41.3 42.5 42.9 44.6 45.6 46.9 46.9 47.6]; % ?

% Pressures
P_H2 =  [2.9 2.9 3.1 3.3 3.30 3.20 3.00 3.0 3.1];                   % psi (gauge)
P_air = [0.2 0.3 0.6 0.7 1.15 1.25 1.35 1.3 1.5];

% CALCULATIONS
% ------------
% Excess Air Coefficient, Lambda
AF = mdot_air./mdot_fuel;
AFs = [];                   % TODO: what is this value?
lambda = Af./AFs;

% Power 
p_in = [];                  % a.k.a. "load" (power delivered to resistor bank)
p_stack = [];
p_access = [];

figure(1)
plot(load,i_load,load,i_stack);
title('Current as a Function of Load');
xlabel('Load []'); ylabel('Current []');
legend('I_{load}','I_{stack}'); plotfixer(); grid on

figure(2) 
plot(load,v_load,load,v_stack);
title('Potential as a Function of Load');
xlabel('Load []'); ylabel('Potential []');
legend('V_{load}','V_{stack}'); plotfixer();grid on

figure(3)
plot(load,p_stack,load,p_access);
title('Stack and Accessory Power as a Function of Load');
xlabel('Load []'); ylabel('Power []');
legend('P_{stack}','P_{accessory}'); plotfixer();grid on

figure(4)
plot(load,mdot_H2,load,mdot_air);
title('Mass Flow Rate as a Function of Load');
xlabel('Load []'); ylabel('Mass Flow Rate []');
legend('mdot_{H}','mdot_{air}'); plotfixer();grid on


% ---------------------------
% Part 2: Reduced-Data Plots
% ---------------------------
% NOTE: Account for stack and load data so that system losses do not effect
% calculation but stack losses do. 

% % STEP 1: Calculate Gibbs Free Energy of Reaction, dG
% % TODO: update this code to work with appropriate T,P,lambda from above
% for Ti = 1:length(T)
%     for pi = 1:length(Ptotal)
%         [etaPres(pi,Ti), pctVapPres(pi,Ti) ,delGPres(pi,Ti),~] ...
%             = PEMstoich(lambda,T(Ti),Ptotal(pi));
%     end
% end

% STEP 2: First Law Efficiency
eta_firstLaw_LHV = dG/LHV;

% STEP 3: Second Law Efficiency (Exergy Efficiency)
% eta_secondLaw = Wdot_net/ (mdot_fuel * deltaG0T) 
% where deltaG0T = delta G of rxn at P0 = 1 bar and T
eta_secondLaw = Wdot_net/(mdot_H2 * dG);
p_loss = [];            % rate of irreversibility


