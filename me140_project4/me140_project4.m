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

N_TO_O = 3.7619;        % Engineering Air Molar Mass Ratio of Nitrogen to Oxygen

% Molar Masses
MM_c = 12;
MM_h = 1;
MM_o = 16;
MM_n = 14;
MM_h2o = 2*MM_h + MM_o;
MM_air = 28.85;

% ----------------------------------------------------
% Part 1 & 2: Efficiency of PEM Fuel Cells Found 3 Ways
% ----------- then, varying lambda & presssure
% -----------------------------------------------------
% ASSUME: isothermal, isobaric i.e. reversible
% USE: First- Law Effiency, eta = (-m_reactants*dg_react)/(mfuel*HV) where HV = LHV or HHV
% SOURCE: LEC 8, SLIDE 13
npts = 100;
HHV_h2 = 141.8*10^6;                    % J/kg,  Higher Heating Value   
LHV_h2 = 120.0*10^6;                    % J/kg,  Lower Heating Value  
% ------------------------------------------
% UNCOMMENT FOR PART 1:
T = linspace(25+C_TO_K,1000+C_TO_K,npts);
lambda = 2;                             % Equivalence Ratio(ASSUME: 100% excess air)        
Patm = 101.3*KPA_TO_PA;                 % Pa,     Preact = Pprod = Patm 
% ------------------------------------------

iterations =0;
eta = zeros(size(T));
pctVap = zeros(size(T));
delG = zeros(size(T));
for i = 1:length(T) %loop temperature (cols)
    [eta(i), pctVap(i),delG(i)] = PEMstoich(lambda,T(i),Patm);
    iterations = iterations + 1;
end
iterations
%PEMstoich assumes per mol of h2, 1mol h2 burned
mass_h2 = 1*(MM_h*2)*G_TO_KG;
eta_HHV = -delG / (HHV_h2 * mass_h2);
eta_LHV = -delG / (LHV_h2 * mass_h2);

eta_carnot = carnotEff(T,T(1));      % ASSUME: Tcold = 25 degrees C

figure(1);
plot(T,eta_HHV,'b--', T,eta_LHV,'m--',T,eta,'g-', T,eta_carnot,'c');
legend('\eta_{HHV}','\eta_{LHV}','\eta_{Mixed Liquid and Gas}','\eta_{Carnot}', 'Location', 'Best');
xlabel('Temperature [K]');
ylabel('Maximum 1st Law Efficiency');
plotfixer();
grid on

% % UNCOMMENT FOR PART 2a (varying lambda)
T_C = [80 220 650 800];
T = T_C + C_TO_K;
lambda = linspace(1,10,npts);       % (Comment back in for Part 2)         
Patm = 101.3*KPA_TO_PA;             % Pa,     Preact = Pprod = Patm 

for Ti = 1:length(T)
    for li = 1:length(lambda)
        [etaLambda(li,Ti), pctVapLambda(li,Ti) ,delGLambda(li,Ti)] ...
            = PEMstoich(lambda(li),T(Ti),Patm);
    end
end
mass_h2 = 1* (MM_h*2)*G_TO_KG;
delH_LHV = LHV_h2 * mass_h2;
etaLambda_LHV = -delGLambda/delH_LHV;
figure(2);
plot(lambda,etaLambda_LHV);
legend('80C','220C','650C','800C','Location','Best');
xlabel('Excess air coefficient \lambda');
ylabel('Efficiency on LHV basis \eta');
plotfixer();
grid on

%%part2.1 plot%%

% % UNCOMMENT FOR PART 2b (varying Patm)
T_C = [80 220 650 800];
T = T_C + C_TO_K;
lambda = 2;                         % Equivalence Ratio(ASSUME: 100% excess air)     
Patm = linspace(101.3*KPA_TO_PA,4053*KPA_TO_PA,npts); % Pa, (Comment back in for Part 2)

for Ti = 1:length(T)
    for pi = 1:length(Patm)
        [etaPres(pi,Ti), pctVapPres(pi,Ti) ,delGPres(pi,Ti)] ...
            = PEMstoich(lambda,T(Ti),Patm(pi));
    end
end

etaPres_LHV = -delGPres/delH_LHV;
figure(3);
plot(Patm/101300,etaPres_LHV);
legend('80C','220C','650C','800C','Location','Best');
xlabel('Pressure - Bar');
ylabel('Efficiency on LHV basis \eta');
plotfixer();
grid on
return

% figure(2)
% plot(T,mol_h2ovap);

% figure(3)
% plot(T,dh,'b',T,hreact,'g');
% plotfixer();

%% Part 3
% what humidity necesary for inlet air to obtain saturated exit?
% below certain temp, condensate forms, so add no water.
% plot inlet air humidity vs T 25-100C

% questions:
% must we take into account the diffusion thru membrane? -> don't need to
% worry about gas diffusion through MEA membrane
lambda = 2; % as before
Ptotal = Patm;
% find psat at exit based on temp, 
T = linspace(25,100,npts);
psat = PsatW(T + C_TO_K);

% find mole fraction of water in products
y_h2o = psat./Ptotal;
mol_out = (mol_o2_prod + mol_h2o + mol_n2);
mol_h2o_sat = mol_out*y_h2o;

% if less than what is formed, add the difference to dry air reagent
alpha = mol_h2o_sat - mol_h2o;
alpha(alpha<0) = 0;
y_h2o_react = alpha./(mol_o2_react + mol_n2 + alpha);  
Pv_react = Ptotal*y_h2o_react;
Pv_react(Pv_react>psat) = psat(Pv_react>psat); % if Pv > psat, Pv = psat 
hum_rel = Pv_react./psat;

% plot relative humidity
% plot(T,alpha,T,omega2,T,hum_rel); - DELETE
% legend('Moles of H2O to Add','Relative Humidity, outlet?') - DELETE
figure(2);
plot(T,hum_rel)
legend('Relative Humidity of Input Air');
xlabel('Temperature - K');
ylabel('Relative Humidity %');
plotfixer();



% strategy for part 3:
% calculate # mols of water in products at saturation (from mole fraction =
% psat/p)--> mol_h2o_sat
% then do atom balance - add water on reactant side - how many mols of
% reactant water would u need to get result 
% Pv/psat of products should be 1, then can use beta = alpha + 1 to relate
% prods and reactants and find Pv of reactants w/ alpha, then just find
% relative humidity = Pv/psat


% % DELETE - FROM BEFORE
% % find partial pressure of products
% y_h2o_prod = mol_h2o/(mol_o2_prod + mol_h2o + mol_n2);
% Pv = y_h2o_prod*Ptotal;
% Pv(Pv>psat) = psat(Pv>psat); % if Pv > psat, Pv = psat 
% figure()
% plot(T, Pv./psat);
% title('product rel hum - should be 1');

%hum_rel = Pv./psat; 
%hum_rel(Pv>psat) = 0; 
%omega = Pv./(Ptotal-Pv)*(MM_h2o)/(MM_air); % formula from lecture does not seem to work.
%omega2 = alpha*(MM_h2o)/(mol_o2_rxn*MM_o*2 + mol_n2*MM_n*2); 

