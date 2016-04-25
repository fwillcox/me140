% ME140 PROJECT 4: FUEL CELLS
% ----------------------------
% FILENAME: me140_project4.m
% Jon Renslo, Emily Bohl, Frankie Willcox, Natasha Berk, Kendall Fagan
% 4/15/16 - Created Jon Renslo

close all; clear; clc;

% Constants
G_TO_KG = 10^-3;
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
% ASSUME: isothermal, isobaric i.e. reversible
% USE: First- Law Effiency, eta = (-m_reactants*dg_rxn)/(mfuel*HV) where HV = LHV or HHV
% SOURCE: LEC 8, SLIDE 13

npts = 100;
T = linspace(25+C_TO_K,1000+C_TO_K,npts);
lambda = 2;                         % Equivalence Ratio(ASSUME: 100% excess air)
HHV_h2 = 141800;                    % kj/kg,  Higher Heating Value   
LHV_h2 = 120000;                    % kj/kg,  Lower Heating Value           
Patm = 101.3*KPA_TO_PA;             % Pa,     Preact = Pprod = Patm

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
<<<<<<< HEAD
gprod_LHV = gEng(T,Patm,'h2ovap',mol_h2o) + gEng(T,Patm,'o2',mol_o2_prod) + gEng(T,Patm,'n2',mol_n2); % J, Gibbs Free Energy 
gprod_HHV = gEng(T,Patm,'h2o',mol_h2o) + gEng(T,Patm,'o2',mol_o2_prod) + gEng(T,Patm,'n2',mol_n2);    
greact = gEng(T,Patm,'h2',mol_h2) + gEng(T,Patm,'o2',mol_o2_rxn) + gEng(T,Patm,'n2',mol_n2);
=======
gprod_LHV = gEng(T,P,'h2ovap',mol_h2o) + gEng(T,P,'o2',mol_o2_prod) + gEng(T,P,'n2',mol_n2); % J, Gibbs Free Energy 
gprod_HHV = gEng(T,P,'h2o',mol_h2o) + gEng(T,P,'o2',mol_o2_prod) + gEng(T,P,'n2',mol_n2);    
greact = gEng(T,P,'h2',mol_h2) + gEng(T,P,'o2',mol_o2_rxn) + gEng(T,P,'n2',mol_n2);
>>>>>>> master

% Account for Gas/Liquid Mixture
% SOURCE: LEC 8 Slide 24, LEC 9, Slide 29
% APPROACH: (1) Assume beta=1, let Pv=Psat (2) Solve for Ptotal
% -------- If Pv < Psat, all vapor. If Pv > Psat, must be some liquid.

beta = 1;               % ASSUME: all vapor
Ptotal = Patm;
Psat = PsatW(T);
Pv = Ptotal*(beta./(beta + 0.5.*(gamma(T)-1) +0.5.*gamma(T).*N_TO_O ));

for i = 1:length(Psat)
    if Pv(i) < Psat(i)
        % All H2O is vapor (beta = 1)
        beta = 1;
        Pv(i) = Ptotal*( beta / ( beta + 0.5*((mol_h2o-beta)-1) +0.5*(mol_h2o-beta)*N_TO_O ) );
    else
        % Some H2O is vapor, some liquid (beta not = 1)
        % LET: Pv = Psat, solve for beta
        Pv(i) = PsatW(i);
        gamma = mol_h2o - beta;
        syms b
        solve(Pv(i)/Ptotal == b/(b + 0.5*((mol_h2o-b)-1) +0.5*(mol_h2o-b)*N_TO_O ) ,b);
        
    end
end
%comment

delG_HHV = gprod_HHV - greact;
delG_LHV = gprod_LHV - greact;
% delG_liquidGasMix = ???
eta_HHV = -delG_HHV / (HHV_h2 * mass_h2 * KJ_TO_J);
eta_LHV = -delG_LHV / (LHV_h2 * mass_h2 * KJ_TO_J);




% Y_h2o = mol_h2o/mol_total;
% eta_liquidGasMix = ??? todo find g actual based on liquid water mixture
% eta_carnot = 1- Tlow/Thigh;

figure();
plot(T,eta_HHV,T,eta_LHV);
legend('\eta_{HHV}','\eta_{LHV}');
xlabel('Temperature K');
ylabel('Maximum 1st Law Efficiency');
plotfixer();


% Part 2

% Plot eta_LHV vs lambda (1:10), (p = 1atm)
%      eta_LHV vs pressure (1:40atm), (lambda = 2)
%       each at 80 220 650 800C

%good luck
% more luck
% even more


