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
MM_air = 28.85;

% ----------------------------------------------------
% Part 1 & 2: Efficiency of PEM Fuel Cells Found 3 Ways
% ----------- then, varrying lambda & presssure
% -----------------------------------------------------
% ASSUME: isothermal, isobaric i.e. reversible
% USE: First- Law Effiency, eta = (-m_reactants*dg_rxn)/(mfuel*HV) where HV = LHV or HHV
% SOURCE: LEC 8, SLIDE 13
npts = 100;
HHV_h2 = 141800;                    % kj/kg,  Higher Heating Value   
LHV_h2 = 120000;                    % kj/kg,  Lower Heating Value  

% ------------------------------------------
% UNCOMMENT FOR PART 1:
T = linspace(25+C_TO_K,1000+C_TO_K,npts);
lambda = 2;                         % Equivalence Ratio(ASSUME: 100% excess air)        
Patm = 101.3*KPA_TO_PA;             % Pa,     Preact = Pprod = Patm 

% % UNCOMMENT FOR PART 2a (varrying lambda)
% T_C = [80 220 650 800];
% T = T_C + C_TO_K;
% lambda = linspace(1,10,npts);       % (Comment back in for Part 2)         
% Patm = 101.3*KPA_TO_PA;             % Pa,     Preact = Pprod = Patm 

% % UNCOMMENT FOR PART 2b (varrying Patm)
% T_C = [80 220 650 800];
% T = T_C + C_TO_K;
% lambda = 2;                         % Equivalence Ratio(ASSUME: 100% excess air)     
% Patm = linspace(101.3*KPA_TO_PA,4053*KPA_TO_PA,npts); % Pa, (Comment back in for Part 2)
% ------------------------------------------

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
gprod_LHV = gEng(T,Patm,'h2ovap',mol_h2o) + gEng(T,Patm,'o2',mol_o2_prod) + gEng(T,Patm,'n2',mol_n2); % J, Gibbs Free Energy 
gprod_HHV = gEng(T,Patm,'h2o',mol_h2o) + gEng(T,Patm,'o2',mol_o2_prod) + gEng(T,Patm,'n2',mol_n2);    
greact = gEng(T,Patm,'h2',mol_h2) + gEng(T,Patm,'o2',mol_o2_rxn) + gEng(T,Patm,'n2',mol_n2);

% Account for Gas/Liquid Mixture
% SOURCE: LEC 8 Slide 24, LEC 9, Slide 29
% APPROACH: (1) Assume beta=1, let Pv=Psat (2) Solve for Ptotal
% -------- If Pv < Psat, all vapor. If Pv > Psat, must be some liquid.
beta = 1;               % ASSUME: all vapor
Ptotal = Patm;
Psat = PsatW(T);
gamma = mol_h2o - beta;
Pv_guess = Ptotal*(beta./(beta + 0.5.*(gamma-1) +0.5.*gamma.*N_TO_O ));
Pv_h2o = Psat;
eta_carnot = carnotEff(T,T(1));      % ASSUME: Tcold = 25 degrees C

iterations =0;
mol_h2ovap=zeros(size(Psat));
for i = 1:length(Psat)
    if Pv_guess < Psat(i)
        % All H2O is vapor (beta = 1)
        mol_h2ovap(i) = beta;
        mol_h2oliq(i) = 0;
        Pv_h2o(i) = Pv_guess;
        
    else
        % Some H2O is vapor, some liquid (beta not = 1)
        % LET: Pv = Psat, solve for beta
        Pv_h2o(i) = Psat(i);
        mol_h2ovap(i) = (mol_o2_prod + mol_n2)*Pv_h2o(i)/(Ptotal-Pv_h2o(i)); % beta
        mol_h2oliq(i) = mol_h2o - mol_h2ovap(i);
        i
    end

% DOUBLE CHECK THIS
mol_total = mol_o2_prod + mol_n2 + mol_h2o;
y_h2ovap = mol_h2o/mol_total;
y_o2 = mol_o2_prod/mol_total;
y_n2 = mol_n2/mol_total;

gprod_LHV_mix(i) = ...
      gEng(T(i), Patm*y_h2ovap, 'h2ovap', mol_h2ovap(i))...
    + gEng(T(i), Patm,          'h2o',    mol_h2oliq(i))...
    + gEng(T(i), Patm*y_o2,     'o2',     mol_o2_prod)...   
    + gEng(T(i), Patm*y_n2,     'n2',     mol_n2);

delG_mix(i) = gprod_LHV_mix(i) - greact(i);    
hprod(i) = ...
      hEng(T(i),'h2ovap', mol_h2ovap(i))...
    + hEng(T(i),'h2o',    mol_h2oliq(i))...
    + hEng(T(i),'o2',     mol_o2_prod)...
    + hEng(T(i),'n2',     mol_n2);
hreact(i) = ...
      hEng(T(i),'h2',     mol_h2)... 
    + hEng(T(i),'o2',     mol_o2_rxn)...
    + hEng(T(i),'n2',     mol_n2);
dh(i) = hprod(i) - hreact(i);
eta_mix(i) = -delG_mix(i)/ dh(i);
iterations = iterations + 1;
end
iterations

delG_HHV = gprod_HHV - greact;
delG_LHV = gprod_LHV - greact;
eta_HHV = -delG_HHV / (HHV_h2 * mass_h2 * KJ_TO_J);
eta_LHV = -delG_LHV / (LHV_h2 * mass_h2 * KJ_TO_J);

figure(1);
plot(T,eta_HHV,'b--', T,eta_LHV,'m--', T,eta_mix,'go', T,eta_carnot,'c');
legend('\eta_{HHV}','\eta_{LHV}','\eta_{Mixed Liquid and Gas}','\eta_{Carnot}');
xlabel('Temperature [K]');
ylabel('Maximum 1st Law Efficiency');
plotfixer();


% %% Part 3
% % what humidity necesarry in inlet air to obtain saturated exit?
% % below certain temp, condensate forms, so add no water.
% % plot inlet air humidity vs T 25-100C
% 
% % questions:
% % must we take into account the diffusion thru membrane?
% lambda = 2; %as before
% Ptotal = Patm;
% % find psat at exit based on temp, 
% T = linspace(25,100,npts);
% psat = PsatW(T+273);
% % find mole fraction of water
% y_h2o = psat./Ptotal;
% y_h2o_prod = mol_h2o/(mol_o2_prod + mol_h2o + mol_n2);
% mol_out = (mol_o2_prod + mol_h2o + mol_n2);
% mol_h2o_sat = mol_out*y_h2o;
% Pv = y_h2o_prod*Ptotal;
% Pv(psat>Pv) = psat(psat>Pv);
% 
% % if less than what is formed, add the difference to dry air reagent
% omega = Pv./(Ptotal-Pv)*(MM_h2o)/(MM_air); %formula from lecture does not seem to work.
% diff = mol_h2o_sat - mol_h2o;
% diff(diff<0) = 0;
% omega2 = diff*(MM_h2o)/(mol_o2_rxn*MM_o*2 + mol_n2*MM_n*2);
% %convert mol fraction to humidity
% plot(T,diff,T,omega2);
% legend('Moles of H2O to Add','Absolute Humidity')


    


