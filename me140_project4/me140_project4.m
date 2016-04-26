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

% % UNCOMMENT FOR PART 2a (varying lambda)
% T_C = [80 220 650 800];
% T = T_C + C_TO_K;
% lambda = linspace(1,10,npts);       % (Comment back in for Part 2)         
% Patm = 101.3*KPA_TO_PA;             % Pa,     Preact = Pprod = Patm 

% % UNCOMMENT FOR PART 2b (varying Patm)
% T_C = [80 220 650 800];
% T = T_C + C_TO_K;
% lambda = 2;                         % Equivalence Ratio(ASSUME: 100% excess air)     
% Patm = linspace(101.3*KPA_TO_PA,4053*KPA_TO_PA,npts); % Pa, (Comment back in for Part 2)
% ------------------------------------------

% assume 1 mol of h2 combusted --> 4.76/2 mol air
% 139 g total
mol_h2 = 1;
mass_h2 = mol_h2*(MM_h*2)*G_TO_KG;
mol_air = 4.76*lambda/2*mol_h2;
mol_o2_react = mol_air/4.76;

mol_n2 = mol_air*3.76/4.76;
mol_o2_prod = 0.5*(lambda-mol_h2)*mol_o2_react;  

mol_h2o = mol_h2;

% Check Mass Balance
% mass_react = mass_h2 + mol_air*MM_air*G_TO_KG;
% mass_prod = (mol_o2_prod*2*MM_o*G_TO_KG) + (mol_n2*2*MM_n*G_TO_KG  ... 
% +mol_h2o*MM_h2o*G_TO_KG);

% Account for Gas/Liquid Mixture
% SOURCE: LEC 8 Slide 24, LEC 9, Slide 29
% APPROACH: (1) Assume beta=1, let Pv=Psat (2) Solve for Ptotal
% -------- If Pv < Psat, all vapor. If Pv > Psat, must be some liquid.
beta = mol_h2;               % ASSUME: all vapor
Ptotal = Patm;
Psat = PsatW(T);
Pv_guess = Ptotal*(beta./(beta + 0.5.*(lambda-1) +0.5.*lambda.*N_TO_O ));
Pv_h2o = zeros(size(Psat));
mol_h2oliq = zeros(size(Psat));
iterations =0;
mol_h2ovap = zeros(size(Psat));
mol_total_prod = zeros(size(Psat));
y_h2ovap = zeros(size(Psat));
y_o2_prod = zeros(size(Psat));
y_n2_prod = zeros(size(Psat));
greact = zeros(size(Psat));
gprod = zeros(size(Psat));
delG = zeros(size(Psat));
hprod = zeros(size(Psat));
hreact = zeros(size(Psat));
dh = zeros(size(Psat));
eta_mix = zeros(size(Psat));

mol_total_react = mol_o2_react + mol_n2;
y_o2_react = mol_o2_react /mol_total_react;
y_n2_react = mol_n2       /mol_total_react;
y_h2_react = mol_h2       /mol_h2; 

for i = 1:length(Psat)
    if Pv_guess < Psat(i)
        % All H2O is vapor (beta = 1)
        mol_h2ovap(i) = beta;
        mol_h2oliq(i) = 0;
        Pv_h2o(i) = Pv_guess;
    else % i = 1-10
        % Some H2O is vapor, some liquid (beta not = 1)
        % LET: Pv = Psat, solve for beta
        Pv_h2o(i) = Psat(i);
        mol_h2ovap(i) = (mol_o2_prod + mol_n2)*Pv_h2o(i)/(Ptotal-Pv_h2o(i)); % beta
        mol_h2oliq(i) = mol_h2o - mol_h2ovap(i);
    end

    % because membrane separates h2 from air, partial pressures are
    % separate
% 
%     greact(i) = gEng(T(i),Patm*y_h2_react,'h2',mol_h2) ...
%               + gEng(T(i),Patm*y_o2_react,'o2',mol_o2_react) ...
%               + gEng(T(i),Patm*y_n2_react,'n2',mol_n2);
% 
%     gprod(i) = ...
%           gEng(T(i), Patm*y_h2ovap(i), 'h2ovap', mol_h2ovap(i))...
%         + gEng(T(i), Patm,          'h2o',       mol_h2oliq(i))...
%         + gEng(T(i), Patm*y_o2_prod(i),     'o2',     mol_o2_prod)...   
%         + gEng(T(i), Patm*y_n2_prod(i),     'n2',     mol_n2);
% 
%     delG(i) = gprod(i) - greact(i);    
%     hprod(i) = ...
%           hEng(T(i),'h2ovap', mol_h2ovap(i))...
%         + hEng(T(i),'h2o',    mol_h2oliq(i))...
%         + hEng(T(i),'o2',     mol_o2_prod)...
%         + hEng(T(i),'n2',     mol_n2);
%     hreact(i) = ...
%           hEng(T(i),'h2',     mol_h2)... 
%         + hEng(T(i),'o2',     mol_o2_react)...
%         + hEng(T(i),'n2',     mol_n2);
%     dh(i) = hprod(i) - hreact(i);
% 
%     eta_mix(i) = delG(i)/ dh(i);

    iterations = iterations + 1;
end
iterations

mol_total_prod  = mol_o2_prod + mol_n2 + mol_h2ovap;
y_h2ovap = mol_h2ovap./mol_total_prod;
y_o2_prod = mol_o2_prod./mol_total_prod;
y_n2_prod = mol_n2./mol_total_prod;

% DOUBLE CHECK THIS
mol_total = mol_o2_prod + mol_n2 + mol_h2ovap;
y_h2ovap = mol_h2ovap ./ mol_total;
y_o2 = mol_o2_prod ./ mol_total;
y_n2 = mol_n2 ./ mol_total;

greact = gEng(T,Patm,'h2',mol_h2) ...
    + gEng(T,Patm .* y_o2,'o2',mol_o2_react) ...
    + gEng(T,Patm .* y_n2,'n2',mol_n2);

gprod = ...
      gEng(T, Patm.*y_h2ovap,   'h2ovap', mol_h2ovap)...
    + gEng(T, Patm,             'h2o',    mol_h2oliq)...
    + gEng(T, Patm.*y_o2,       'o2',     mol_o2_prod)...   
    + gEng(T, Patm.*y_n2,       'n2',     mol_n2);

delG = gprod - greact;    
hprod = ...
      hEng(T,'h2ovap', mol_h2ovap)...
    + hEng(T,'h2o',    mol_h2oliq)...
    + hEng(T,'o2',     mol_o2_prod)...
    + hEng(T,'n2',     mol_n2);
hreact = ...
      hEng(T,'h2',     mol_h2)... 
    + hEng(T,'o2',     mol_o2_react)...
    + hEng(T,'n2',     mol_n2);
dh = hprod - hreact;

eta_mix = delG ./ dh;



eta_HHV = -delG / (HHV_h2 * mass_h2);
eta_LHV = -delG / (LHV_h2 * mass_h2);

eta_carnot = carnotEff(T,T(1));      % ASSUME: Tcold = 25 degrees C


figure(1);
plot(T,eta_HHV,'b--', T,eta_LHV,'m--',T,eta_mix,'g-', T,eta_carnot,'c');
legend('\eta_{HHV}','\eta_{LHV}','\eta_{Mixed Liquid and Gas}','\eta_{Carnot}', 'Location', 'Best');
xlabel('Temperature [K]');
ylabel('Maximum 1st Law Efficiency');
plotfixer();

% figure(2)
% plot(T,mol_h2ovap);

% figure(3)
% plot(T,dh,'b',T,hreact,'g');
% plotfixer();




%% Part 3
% what humidity necesarry in inlet air to obtain saturated exit?
% below certain temp, condensate forms, so add no water.
% plot inlet air humidity vs T 25-100C

% questions:
% must we take into account the diffusion thru membrane? -> don't need to
% worry about gas diffusion through MEA membrane
lambda = 2; %as before
Ptotal = Patm;
% find psat at exit based on temp, 
T = linspace(25,100,npts);
psat = PsatW(T+273);

% find mole fraction of water in products
y_h2o = psat./Ptotal;
mol_out = (mol_o2_prod + mol_h2o + mol_n2);
mol_h2o_sat = mol_out*y_h2o;

y_h2o_prod = mol_h2o/(mol_o2_prod + mol_h2o + mol_n2);


% if less than what is formed, add the difference to dry air reagent
% omega = Pv./(Ptotal-Pv)*(MM_h2o)/(MM_air); % DELETE? - formula from lecture does not seem to work.
alpha = mol_h2o_sat - mol_h2o;
alpha(alpha<0) = 0;
y_h2o_react = alpha./(mol_o2_rxn + mol_n2 + alpha);  
Pv_react = Ptotal*y_h2o_react;
Pv_react(Pv_react>psat) = psat(Pv_react>psat); % if Pv > psat, Pv = psat 
hum_rel = Pv_react./psat;
hum_rel(Pv_react>psat) = 0; 

% omega2 = alpha*(MM_h2o)/(mol_o2_rxn*MM_o*2 + mol_n2*MM_n*2); - DELETE
%convert mol fraction to humidity
% plot(T,alpha,T,omega2,T,hum_rel);
% legend('Moles of H2O to Add','Relative Humidity, outlet?')
plot(T,hum_rel)
legend('Relative Humidity of Input Air');
xlabel('Temperature - K');
ylabel('Relative Humidity %');
plotfixer();

% calculate # mols of water in products at saturation (from mole fraction =
% psat/p)--> mol_h2o_sat
% then do atom balance - add water on reactant side - how many mols of
% reactant water would u need to get result 
% Pv/psat of products should be 1, then can use beta = alpha + 1 to relate
% prods and reactants and find Pv of reactants w/ alpha, then just find
% Pv/psat
% relative humidity is just Pv/psat



% DELETE
% % find partial pressure of products
% Pv = y_h2o_prod*Ptotal;
% Pv(Pv>psat) = psat(Pv>psat); % if Pv > psat, Pv = psat 

%hum_rel = Pv./psat; % DELETE
%hum_rel(Pv>psat) = 0; % DELETE

% DELETE
% figure();
% plot(T,Pv./psat);
% title('Pv/psat');
    

