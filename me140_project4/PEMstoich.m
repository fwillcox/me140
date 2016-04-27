% PEMstoich.m
% 4-26-16 Created Jon Renslo
function [eta, pctVap,delG, specs] = PEMstoich(lambda,T,Ptotal,alpha)
% Calculates the stoichiometry for a PEM fuel cell, dry air and dry h2
% Returns the efficiency and the % of water vapor in the products by mass

% do we want to return a mixture vector also?
if(~exist('alpha','var')) 
    alpha = 0; 
end

% all return values per mol of fuel burned (assuming 1 mol here)
specs = Spec(); %class initialization
mol_h2 = 1; 
mol_air = 4.76*lambda/2*mol_h2;
mol_o2_react = mol_air/4.76;
N_TO_O = 3.7619;        % Engineering Air Molar Mass Ratio of Nitrogen to Oxygen

mol_n2 = mol_air*3.76/4.76;
mol_o2_prod = 0.5*(lambda-1).*mol_h2;  
%double check o2prod? should be *mol_h2?

mol_h2o = mol_h2 + alpha;

beta = mol_h2o;               % ASSUME: all vapor
% Ptotal = Patm; from before restructuring
Psat = PsatW(T);
Pv_guess = Ptotal*(beta./(beta + 0.5.*(lambda-1) +0.5.*lambda.*N_TO_O ));

mol_total_react = mol_o2_react + mol_n2 + alpha;
y_o2_react = mol_o2_react /mol_total_react;
y_n2_react = mol_n2       /mol_total_react;
y_h2o_react = alpha       /mol_total_react;

% because membrane separates h2 from air, partial pressures are separate

if Pv_guess < Psat
    % All H2O is vapor (beta = 1)
    mol_h2ovap = beta;
    mol_h2oliq = 0;
else % i = 1-10
    % Some H2O is vapor, some liquid (beta not = 1)
    % LET: Pv = Psat, solve for beta
    Pv_h2o = Psat;
    mol_h2ovap = (mol_o2_prod + mol_n2)*Pv_h2o/(Ptotal-Pv_h2o); %  = beta
    mol_h2oliq = mol_h2o - mol_h2ovap;
end

pctVap = mol_h2ovap./(mol_h2o);

% With the moles of liquid and gas water products, calculate mole fractions
% and deltaG and deltaH
mol_total_prod  = mol_o2_prod + mol_n2 + mol_h2ovap;
y_h2ovap  = mol_h2ovap  ./ mol_total_prod;
y_o2_prod = mol_o2_prod ./ mol_total_prod;
y_n2_prod = mol_n2      ./ mol_total_prod;

greact = gEng(T,Ptotal,'h2',mol_h2) ...
    + gEng(T,Ptotal .* y_o2_react,'o2',mol_o2_react) ...
    + gEng(T,Ptotal .* y_n2_react,'n2',mol_n2);
if(alpha ~=0)
    greact = greact + gEng(T,Ptotal .* y_h2o_react,'h2ovap',alpha);
end

gprod = ...
      gEng(T, Ptotal.*y_h2ovap,   'h2ovap', mol_h2ovap)...
    + gEng(T, Ptotal,             'h2o',    mol_h2oliq)...
    + gEng(T, Ptotal.*y_o2_prod,       'o2',     mol_o2_prod)...   
    + gEng(T, Ptotal.*y_n2_prod,       'n2',     mol_n2);

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

eta = delG ./ dh;
specs.mol_air =         mol_air;
specs.mol_o2_react =         mol_o2_react;
specs.mol_n2 =         mol_n2;
% TODO update Spec to accomodate inlet water? (PLEASE CHECK KENDALL)

end