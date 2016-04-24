% me140_project4.m
% 4-22-16 - Created Jon Renslo
% Script for Project 4 in ME140. Fuel Cells. 
close all; clear; clc;

% Constants
G_TO_KG = -10^3;
KPA_TO_PA = 10^3;
KJ_TO_J = 10^3;
C_TO_K = 273.15;

N_TO_O = 79/21;         % Engineering Air Molar Mass Ratio of Nitrogen to Oxygen
AIR_TO_H = 4.76;        % Ratio between H combusted and air

% Molar Masses
MM_c = 12;
MM_h = 1;
MM_o = 16;
MM_n = 14;
MM_h2o = 2*MM_h + MM_o;

% -------------------------------------------------
% Part 1: Efficiency of PEM Fuel Cells Found 3 Ways
% -------------------------------------------------
% (1) LHV, (2) HHV, (3) Accounting for liquid/gas mixture

npts = 100;
T = linspace(25+C_TO_K,1000+C_TO_K,npts);
lambda = 2;                         % Equivalence Ratio(ASSUME: 100% excess air)
HHV_h2 = 141800;                    % kj/kg,  Higher Heating Value   
LHV_h2 = 120000;                    % kj/kg,  Lower Heating Value           
P = 101.3*KPA_TO_PA;                % Pa,     Preact = Pprod = Patm

mol_h2 = 1;                         % (ASSUME: 1 mol H2-->4.76/2 mol air = 139 g)
mass_h2 = mol_h2 * 2*MM_h * G_TO_KG; 
mol_air = AIR_TO_H * lambda / 2;
mol_o2_rxn = mol_air / AIR_TO_H;
mol_n2 = mol_air * N_TO_O / AIR_TO_H;

mol_h2o = mol_h2;
mol_o2_prod = 0.5*(lambda - mol_h2) * mol_o2_rxn;

% Check Mass Balance
% MM_air = 28.85;
% mass_react = mass_h2 + mol_air*MM_air*G_TO_KG
% mass_prod = (mol_o2_prod*2*MM_o*G_TO_KG) + (mol_n2*2*MM_n*G_TO_KG + ...
% ...mol_h2o*MM_h2o*G_TO_KG);

% Calculate Change in Gibbs Free Energy 
gprod_LHV = gEng(T,P,'h2o_vap',mol_h2o) + gEng(T,P,'o2',mol_o2_prod) + gEng(T,P,'n2',mol_n2); % J, Gibbs Free Energy 
gprod_HHV = gEng(T,P,'h2o',mol_h2o) + gEng(T,P,'o2',mol_o2_prod) + gEng(T,P,'n2',mol_n2);    
greact = gEng(T,P,'h2',mol_h2) + gEng(T,P,'o2',mol_o2_rxn) + gEng(T,P,'n2',mol_n2);

delG_HHV = gprod_HHV - greact;
delG_LHV = gprod_LHV - greact;
% delG_liquidGasMix = ???
eta_HHV = -delG_HHV / (HHV_h2 * mass_h2 * KJ_TO_J);
eta_LHV = -delG_LHV / (LHV_h2 * mass_h2 * KJ_TO_J);
% eta_liquidGasMix = ??? todo find g actual based on liquid water mixture
% eta_ carnot

figure();
plot(T,eta_HHV,T,eta_LHV);
legend('\eta_{HHV}','\eta_{LHV}');
xlabel('Temperature K');
ylabel('Maximum 1st Law Efficiency');
plotfixer();


