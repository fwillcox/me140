% ME 140 Project #5
% FUEL CELL EVALUATION & HYRDOGEN PRODUCTION ANALYSIS
% Frankie Willcox, Jon Renslo, Kendall Fagan, Emily Bohl, Natasha Berke

% Constants

% --------------------------------
% Part 1: Raw Data Plots vs. Load
% --------------------------------
% Load (Power delivered only to resistor bank)
load = [];

% Currents (load & stack)
i_load = [];
i_stack = [];

% Potentials (load & stack)
v_load = [];
v_stack = [];

% Power (load & accessory)
p_stack = [];
p_access = [];

% Mass Flow Rates (H & air)
mdot_H = [];
mdot_air = [];

f1 = figure(1);
plot(load,i_load,load,i_stack);
title('Current as a Function of Load');
xlabel('Load []'); ylabel('Current []');
legend('I_{load}','I_{stack}'); plotfixer(); grid on;
saveas(f1,'1-CurrentbyLoad','png');

f2 = figure(2) 
plot(load,v_load,load,v_stack);
title('Potential as a Function of Load');
xlabel('Load []'); ylabel('Potential []');
legend('V_{load}','V_{stack}'); plotfixer();grid on;
saveas(f2,'2-VbyLoad','png');


f3 = figure(3)
plot(load,p_stack,load,p_access);
title('Stack and Accessory Power as a Function of Load');
xlabel('Load []'); ylabel('Power []');
legend('P_{stack}','P_{accessory}'); plotfixer();grid on;
saveas(f3,'3-PowerbyLoad','png');

f4 = figure(4);
plot(load,mdot_H,load,mdot_air);
title('Mass Flow Rate as a Function of Load');
xlabel('Load []'); ylabel('Mass Flow Rate []');
legend('mdot_{H}','mdot_{air}'); plotfixer();grid on
saveas(f4,'4-massbyload','png');



% ---------------------------
% Part 2: Reduced-Data Plots
% ---------------------------
lambda = [];            % Excess Air Coefficient

% NOTE: Account for stack and load data so that system losses do not effect
% calculation but stack losses do. 
eta_firstlaw_LHV = [];
eta_secondlaw = [];
p_loss = [];            % rate of irreversibility


