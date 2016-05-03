% ME 140 Project #5
% FUEL CELL EVALUATION & HYRDOGEN PRODUCTION ANALYSIS
% Frankie Willcox, Jon Renslo, Kendall Fagan, Emily Bohl, Natasha Berke

% ASSUME:
% (i)  mol_H2 = 1
clear; close all;clc;

global PERMIN_TO_PERSEC PERHR_TO_PERSEC G_PER_KG LHV F N_TO_O SCF_TO_MOLS ...
    C_TO_K PSI_TO_PA MM_h MM_h2 MM_o MM_n MM_h2o MM_air PATM HORSEPOWER_TO_W
defineGlobals();
mol_H2 = 1;
savePlots = 1;
                % 1,2,3,4,5,6,7,8,9
supressplots =   [1,1,0,1];         % supresses plots by section
plotNum = 1;

%% Part A, Section 1
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
T4_C =     [48.0 47.2 48.2 48.9 50.4 51.1 51.2 51.1 51.1];         
T5_C =     [40.7 41.3 42.5 42.9 44.6 45.6 46.9 46.9 47.6];       
T4 = T4_C + C_TO_K;                                                % T4, water into stack
T5 = T5_C + C_TO_K;                                                % T5, water into heat exchanger


% Mass Flow Rates (TODO: check what units the mdots should be in)
mdot_total_scf = [0.75 1.10 1.45 1.81 2.55 3.10 3.30 3.25 3.40];   % [scf/min]
mdot_fuel_scf =  [2.50 6.20 10.5 14.3 18.2 22.0 24.6 25.0 26.1];   % [scf/hr] (standard cubic feet/hour)

mdot_total = mdot_total_scf * SCF_TO_MOLS * PERMIN_TO_PERSEC * ( MM_air / G_PER_KG);  % [kg/s]
mdot_fuel = mdot_fuel_scf * SCF_TO_MOLS * PERHR_TO_PERSEC * (MM_h2  / G_PER_KG);      % [kg/s]	
mdot_h2o = 40 /G_PER_KG;                                                              % [kg/s]

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
if(~supressplots(plotNum))
f1 = figure(1);
plot(p_load,i_load,p_load,i_stack);
title('Current as a Function of Load');
xlabel('Load [Watts]'); ylabel('Current [Amps]');
legend('I_{load}','I_{stack}','Location','best'); plotfixer(); grid on;

f2 = figure(2); 
plot(p_load,v_load,p_load,v_stack);
title('Potential as a Function of Load');
xlabel('Load [Watts]'); ylabel('Potential [Volts]');
legend('V_{load}','V_{stack}','Location','best'); plotfixer();grid on;

f3 = figure(3);
plot(p_load,p_stack,p_load,p_access);
title('Stack and Accessory Power as a Function of Load');
xlabel('Load [Watts]'); ylabel('Power [Watts]');
legend('P_{stack}','P_{accessory}','Location','best'); plotfixer();grid on;

f4 = figure(4);
plot(p_load,mdot_fuel,p_load,mdot_total);
title('Mass Flow Rate as a Function of Load');
xlabel('Load [Watts]'); ylabel('Mass Flow Rate [kg/s]');
legend('mdot_{H}','mdot_{air}','Location','best'); plotfixer(); grid on

plotNum = plotNum+1;
end

%% Part A, Section 2

% SOURCE: LEC 8, SLIDES 21 & 22

% 1st & 2nd Law Efficiencies (eta_I & eta_II) & Inefficiencies (Idot)
% Stack
[etaI_stack ,etaII_stack, Idot_stack,lambda_stack] = findEtas(mdot_total, mdot_fuel, Ptotal, Pfuel, T2, p_stack);

% Entire System (Load)
[etaI_load ,etaII_load, Idot_load,lambda_load] =    findEtas(mdot_total, mdot_fuel, Ptotal, Pfuel, T2, p_load);
if(~supressplots(plotNum))
f6 = figure(6);
plot(p_load,lambda_load);
title('Air Equivalent as a Function of Load');
xlabel('Load [Watts]'); ylabel('Lambda');
legend('\lambda','Location','best'); plotfixer();grid on

f5 = figure(5);
plot(p_load,etaI_stack,'c',p_load,etaI_load,'bp--',...
    p_load,etaII_stack,'r',p_load,etaII_load,'gp--');
title('Efficiency as a Function of Load');
xlabel('Load [Watts]'); ylabel('Efficiency');
legend('\eta_{I,stack}','\eta_{I,system}',...
    '\eta_{II,stack}','\eta_{II,system}', 'Location','Best'); plotfixer(); grid on;

f7 = figure(7);
plot(p_load,p_stack,'c',p_load,p_load,'bp--');
title('Power Loss/Inefficiences as a Function of Load');
xlabel('Load [Watts]'); ylabel('Power Loss/Inefficiencies, Idot [Watts]');
legend('Idot_{stack}','Idot_{system}','Location','best'); plotfixer();grid on;
plotNum = plotNum+1;
end

%% Part A, Section 3
% Comparing First Law Efficiencies of PEM Fuel Cell with Diesel & Hybrid Engines

% Typical modern Diesel engine (eta_disel = 42%) (chose diesel truck because it's better than a car and worse than a freight ship)
% Source Efficiency: Slide 3, http://www.sae.org/events/gim/presentations/2011/RolandGravel.pdf
% Source Horsepower: https://cumminsengines.com/isx15-heavy-duty-truck-2013#overview
eta_diesel = 0.42;
Wdot_diesel = 400 * HORSEPOWER_TO_W;  % [W]

% Typical gasoline hybrid engine (eta_hybrid = max of 40%)
% Source Efficiency & Horsepower: Toyota Hybrid Vehicles, http://www.toyota-global.com/innovation/environmental_technology/hybrid/
eta_hybrid = 0.40;
Wdot_hybrid = 121 * HORSEPOWER_TO_W;  % [W]

% Calcuate Heat Removal (Qdot)
for i = 1:length(T4)
Qdot_fuelCell(i) = hEng(T4(i),'h2o') - hEng(T5(i),'h2o');
end
Qdot_fuelCell_max = max(Qdot_fuelCell);

% Theoretical Number of Fuel Cells Needed
num_fuelCells_diesel = Wdot_diesel ./ Qdot_fuelCell_max;
num_fuelCells_hybrid = Wdot_hybrid ./ Qdot_fuelCell_max;

if(~supressplots(plotNum))
% Overall First Law Efficiency of the PEM Fuel Cell = Stack Efficiency
f8 = figure(8);
plot(p_load, etaI_stack, 'c', p_load, eta_diesel, 'bp--', p_load, eta_hybrid, 'gd');
title('Comparing 1st Law Efficiency: PEM Fuel Cell, Diesel, and Gasoline Hybrid');
xlabel('Load [Watts]'); ylabel('Efficiency, eta_{I}');
legend('eta_{I,stack}','eta_{I,Diesel}', 'eta_{I,Hybrid}','Location','best'); plotfixer(); grid on;

plotNum = plotNum+1;
end

%% Part B, Section 1
% Part B, Section 1 - Emily & Kendall
% Calculating Kp Values
% SOURCE Kp Formula: LECTURE 14, SLIDE 4

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
P_ref = 101300; %This is the pressure defined for standard conditions. 
                % Standard conditions are what we need because that is what 
                % the little zero indicates in the equation for g. 
R_u = 8.314; %Universal gas constant

%G_reaction = G_products - G_reactants
g_SMR = (gEng(T_B1, P_ref, 'co',v_CO_SMR) + gEng(T_B1, P_ref, 'h2',v_H2_SMR)) - ...
    (gEng(T_B1, P_ref, 'h2ovap',v_H2O_SMR) + gEng(T_B1, P_ref, 'ch4',v_CH4_SMR));

g_WGS = (gEng(T_B1, P_ref, 'h2',v_H2_WGS) + gEng(T_B1, P_ref, 'co2',v_CO2_WGS)) - ...
    (gEng(T_B1, P_ref, 'h2ovap',v_H2O_WGS) + gEng(T_B1, P_ref, 'co',v_CO_WGS));

%Lecture 13 - Slide 15
kp_SMR = exp(-g_SMR ./ (R_u .* T_B1)); %increases with temp
kp_WGS = exp(-g_WGS ./ (R_u .* T_B1)); %decrease with temp


%functions for convenience
f_kp_SMR = @(T_B1) exp(-((gEng(T_B1, P_ref, 'co',v_CO_SMR) ...
                        + gEng(T_B1, P_ref, 'h2',v_H2_SMR))  ...
                      - (gEng(T_B1, P_ref, 'h2ovap',v_H2O_SMR) ...
                            + gEng(T_B1, P_ref, 'ch4',v_CH4_SMR)))...
                    ./ (R_u.*T_B1));
                
f_kp_WGS = @(T_B1) exp(-((gEng(T_B1, P_ref, 'h2',v_H2_WGS) ...
                            + gEng(T_B1, P_ref, 'co2',v_CO2_WGS)) ...
                         -(gEng(T_B1, P_ref, 'h2ovap',v_H2O_WGS) ...
                            + gEng(T_B1, P_ref, 'co',v_CO_WGS))) ...
                        ./ (R_u.*T_B1));

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
if(~supressplots(4))
f9 = figure(9);
kpIsOne = ones(size(T_B1));
semilogy(T_B1(i_min_SMR:i_max_SMR), kp_SMR(i_min_SMR:i_max_SMR), ...
         T_B1(i_min_WGS:i_max_WGS), kp_WGS(i_min_WGS:i_max_WGS));
xlabel('Temperature [C]')
ylabel('Equilibrium Constant')
legend('SMR', 'WGS')
title('Part B.1: Equilibrium Constant vs. Temperature')
ylim([0.001,1000]);
plotfixer();grid on
end

%% Part B No. 2
% Find the Equilibrium Composition (Mol Fractions) of the Steam Methane 
% Reformation(SMR) Reaction 
% SOURCE Nernst Atom Balance: LECTURE 14, SLIDE 4 (equation in lower right corner)
npts = 20;
syms nco nch4 nh2 nh2o;

temps = linspace(25,1200,npts);
temps = temps + C_TO_K;
pres = [1,10,100];
soln = zeros(length(temps),4,length(pres));
tic
for i = 1:length(temps)
    for j = 1:length(pres)
        p = pres(j);
        t = temps(i);
        
        eqs = [1  == nco   + nch4;...             carbon atom balance
               10 == nh2*2 + nch4*4 + nh2o*2; ... hydrogen atom balance
               3  == nco   + nh2o;...             oxygen atom balance
               nco.*nh2.^3./(nch4.*nh2o).* ...    Nernst atom balance
                       (p ./ (nco + nch4 + nh2 + nh2o)) ...
                     == f_kp_SMR(t)]; 
        % 4 eq, 4 unknown
        [a,b,c,d] = vpasolve(eqs,[nco,nch4,nh2,nh2o],[1,1,1,1]);
        soln(i,:,j) = double(max(real([a,b,c,d])));
        
    end 
end
toc

semilogy(temps,soln(:,1,1),'b',temps,soln(:,1,2),'--b',temps,soln(:,1,3),'.b',...
        temps,soln(:,2,1),'r',temps,soln(:,2,2),'--r',temps,soln(:,2,3),'.r',...
        temps,soln(:,3,1),'k',temps,soln(:,3,2),'--k',temps,soln(:,3,3),'.k',...
        temps,soln(:,4,1),'m',temps,soln(:,4,2),'--m',temps,soln(:,4,3),'.m');


if(savePlots ==1) 
    if(~supressplots(1))
    saveas(f1,'../plots5/1-CurrentbyLoad','png');
    saveas(f2,'../plots5/2-VbyLoad','png'); 
    saveas(f3,'../plots5/3-PowerbyLoad','png'); 
    saveas(f4,'../plots5/4-massbyload','png'); 
    end
    if(~supressplots(2))
    saveas(f5,'../plots5/5-Eff','png'); 
    saveas(f6,'../plots5/6-lambda','png');
    saveas(f7,'../plots5/7-PowerLoss','png');
    end
    if(~supressplots(3))
    saveas(f8,'../plots5/8-CompareToGasoline','png');
    end
    if(~supressplots(4))
    saveas(f9,'../plots5/9-KeqbyT','png');
    end
end

% Part B No. 3
Find the Equilibrium Composition (Mol Fractions) of the Water-Gas Shift 
(WGS) Reaction 
% Equations we'll need:
 eqs = [       1  == nco2   + nco;...          % carbon atom balance
               4  == nco2*2 + nco + nh2o; ...  % hydrogen atom balance
               6  == nh2*2   + nh2o*2;...      % oxygen atom balance
               nco.*nh2o.^3./(nco2.*nh2).* ... % Nernst atom balance 
                  == f_kp_SMR(t)]; 
        % 4 eq, 4 unknown
        [a,b,c,d] = vpasolve(eqs,[nco,nh2o,nco2,nh2],[1,1,1,1]);


%% Part B No. 4
% Plot exit composition (mol fractions) vs. 3 system stations (Reformer,
% Shift Reactor 1, Shift Reactor 2) 
% Note: do this for 2 Different Assumptions: (1) isothermal, (2) adiabatic
% SMR: CH4 + 3*H2O --> CO + 3*H2
% WGS: CO  + 3*H2O --> CO2 + H2

% NAMING CONVENTIONS: 
% Station Location: 1=Reformer, 2 = 1st Shift Reactor, 3 = 2nd Shift
% Assumption:       iso = isothermal, adi = adiabatic 
% Inlet/Exit:       in = inlet, ex = exit

% Inlet Temperatures 
Tin_iso_C = [800 400 250];    % [C]
Tin_adi_C = [800 0 0];        % [C] TODO: solve for Tin_adi_C(2) & (3)
Tin_iso = Tin_iso_C * C_TO_K; % [K]
Tin_adi = Tin_adi_C * C_TO_K; % [K]

% Exit Temperatures
Tex_iso_C = [800 400 250];
Tex_adi_C = [800 0 0];
Tex_iso = Tex_iso_C * C_TO_K;
Tex_adi = Tex_adi_C * C_TO_K;

% Heat Addition for Isothermal Reaction (Qin, ASSUME: isothermal)
Qin_iso = [0 0 0];             % [MJ/(kg of reactants)]

% Percent Methane Burned to Heat Reformer (pct_CH4, ASSUME: adiabatic)
pct_CH4 = [0]; % Note: only applies to Reformer! Not Shift Reactors!

% Efficiency (eta of entire system: reformer & both shift reactors)
% eta = (LHV_h2*mass_h2) / (LHV_ch4*mass_ch4)
eta_iso = [];
eta_adi = [];





