% PEMstoich.m
% 4-26-16 Created Jon Renslo
function [eta, pctVap,delG] = PEMstoich(lambda,T,Ptotal)
% Calculates the stoichiometry for a PEM fuel cell
% Returns the efficiency and the % of water vapor in the products by mass

% do we want to return a mixture vector also?

% all return values per mol of fuel burned (assuming 1 mol here)
mol_h2 = 1; 
mol_air = 4.76*lambda/2*mol_h2;
mol_o2_react = mol_air/4.76;
N_TO_O = 3.7619;        % Engineering Air Molar Mass Ratio of Nitrogen to Oxygen

mol_n2 = mol_air*3.76/4.76;
mol_o2_prod = 0.5*(lambda-1).*mol_h2;  
%double check o2prod? should be *mol_h2?

mol_h2o = mol_h2;

beta = mol_h2;               % ASSUME: all vapor
% Ptotal = Patm; from before restructuring
Psat = PsatW(T);
Pv_guess = Ptotal*(beta./(beta + 0.5.*(lambda-1) +0.5.*lambda.*N_TO_O ));

mol_total_react = mol_o2_react + mol_n2;
y_o2_react = mol_o2_react /mol_total_react;
y_n2_react = mol_n2       /mol_total_react;
% because membrane separates h2 from air, partial pressures are separate

if Pv_guess < Psat
    % All H2O is vapor (beta = 1)
    mol_h2ovap = beta;
    mol_h2oliq = 0;
else % i = 1-10
    % Some H2O is vapor, some liquid (beta not = 1)
    % LET: Pv = Psat, solve for beta
    Pv_h2o = Psat;
    mol_h2ovap = (mol_o2_prod + mol_n2)*Pv_h2o/(Ptotal-Pv_h2o); % beta
    mol_h2oliq = mol_h2o - mol_h2ovap;
end

pctVap = mol_h2ovap/(mol_h2o);

mol_total_prod  = mol_o2_prod + mol_n2 + mol_h2ovap;%% sizes do not agree
y_h2ovap = mol_h2ovap./mol_total_prod;
y_o2_prod = mol_o2_prod./mol_total_prod;
y_n2_prod = mol_n2./mol_total_prod;

greact = gEng(T,Ptotal,'h2',mol_h2) ...
    + gEng(T,Ptotal .* y_o2_react,'o2',mol_o2_react) ...
    + gEng(T,Ptotal .* y_n2_react,'n2',mol_n2);

gprod = ...
      gEng(T, Ptotal.*y_h2ovap,   'h2ovap', mol_h2ovap)...
    + gEng(T, Ptotal,             'h2o',    mol_h2oliq)...
    + gEng(T, Ptotal.*y_o2_prod,       'o2',     mol_o2_prod)...   
    + gEng(T, Ptotal.*y_n2_prod,       'n2',     mol_n2);

if(~isequal(size(greact),size(gprod)))
    mol_o2_prod = repmat(mol_o2_prod',size(T));
    mol_n2 = repmat(mol_t2',size(T));
end

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
end