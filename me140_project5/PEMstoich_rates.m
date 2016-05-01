% PEMstoich_rates.m
% Calculates the stoichiometry for a PEM fuel cell, dry air and dry h2
% and returns the first and second law efficiencies

function [eta_I, eta_II] = PEMstoich_rates(lambda, T, Ptotal, alpha, mdot_h2, mdot_air, mdot_h2o)

if(~exist('alpha','var'))                      
    alpha = 0; 
end

% Constants
N_TO_O = 79/21;                                 % Engineering Air Molar Mass Ratio of Nitrogen to Oxygen
specs = Spec();                                 % class initialization

% Mols of Each Species
% --------------------
% Both (N2 stays constant)
mol_n2 = mol_air*N_TO_O/(1+N_TO_O);

% Reactants
mol_h2 = 1;                                     % ASSUME: 1 mol of air
mol_air = (1+N_TO_O)*lambda/2*mol_h2;
mol_o2_react = mol_air/(1+N_TO_O);
mol_total_react = mol_o2_react + mol_n2 + alpha;

% Products
mol_o2_prod = 0.5*(lambda-1).*mol_h2;           % TODO: double check o2prod. Should it be *mol_h2?
mol_h2o = mol_h2 + alpha;


% Mol Fractions of Each Species
% ------------------------------
y_o2_react = mol_o2_react /mol_total_react;
y_n2_react = mol_n2       /mol_total_react;
y_h2o_react = alpha       /mol_total_react;


% Vapor Fraction of H2O 
% --------------------------
beta = mol_h2o;                                  %(ASSUME: All H2O is vapor in order to calculate Pv_guess)
Psat = PsatW(T);
Pv_guess = Ptotal*(beta./(beta + 0.5.*(lambda-1) +0.5.*lambda.*N_TO_O ));
if Pv_guess < Psat
    % All H2O is vapor (beta = 1)
    mol_h2ovap = beta;
    mol_h2oliq = 0;
else 
    % Some H2O is vapor, some liquid (beta not = 1)
    Pv_h2o = Psat;                      % LET: Pv = Psat to solve for beta
    y_h2o = Pv_h2o./Ptotal; 
    beta = (4.26 .* y_h2o)./ (1 - y_h2o);
    mol_h2ovap = beta; 
    mol_h2oliq = mol_h2o - mol_h2ovap;
end

pctVap = mol_h2ovap./(mol_h2o);         % Percent of H2O that is vapor (by mass)


% Mol Fractions of Products
% --------------------------
mol_total_prod  = mol_o2_prod + mol_n2 + mol_h2ovap;
y_h2ovap  = mol_h2ovap  ./ mol_total_prod;
y_o2_prod = mol_o2_prod ./ mol_total_prod;
y_n2_prod = mol_n2      ./ mol_total_prod;


% Change in Gibbs Free Energy (dG)
% --------------------------------
% Reactants
Gdot_react = gEng(T,Ptotal,'h2',mol_h2)                    * mdot_h2 ...
           + gEng(T,Ptotal .* y_o2_react,'o2',mol_o2_react)* mdot_air ...
           + gEng(T,Ptotal .* y_n2_react,'n2',mol_n2)      * mdot_air;
if(alpha ~= 0 ) 
    Gdot_react = Gdot_react + gEng(T,Ptotal*y_h2o_react,'h2ovap',alpha)* mdot_h2o; 
end

% Products
Gdot_prod = gEng(T, Ptotal*y_h2ovap,   'h2ovap', mol_h2ovap)* mdot_h2o... %% THIS WAS THE PART 2 BUG FROM PROJECT #4!
      + gEng(T, Ptotal,            'h2o',    mol_h2oliq)    * mdot_h2o...
      + gEng(T, Ptotal.*y_o2_prod, 'o2',     mol_o2_prod)   * mdot_air...   
      + gEng(T, Ptotal.*y_n2_prod, 'n2',     mol_n2)        * mdot_air;

dGdot = Gdot_prod - Gdot_react; 


% Change in Enthalpy (dH)
% -----------------------
% Reactants
Hdot_react = hEng(T,'h2',     mol_h2)      * mdot_h2... 
           + hEng(T,'o2',     mol_o2_react)* mdot_air...
           + hEng(T,'n2',     mol_n2)      * mdot_air;
if(alpha ~= 0 )
    Hdot_react = Hdot_react + hEng(T,'h2ovap',alpha) * mdot_h2o; 
end

% Products
Hdot_prod = hEng(T,'h2ovap', mol_h2ovap) * mdot_h2o...
          + hEng(T,'h2o',    mol_h2oliq) * mdot_h2o...
          + hEng(T,'o2',     mol_o2_prod)* mdot_air...
          + hEng(T,'n2',     mol_n2)     * mdot_air;

dHdot = Hdot_prod - Hdot_react;


% Rate of Work (Wdot)
% -------------------
Ndot = mol_H2 *mdot_H2;         % TODO: OK to assume mol_H2 = 1??
Idot  = Ndot * F;              
Wdot  = -dGdot - Idot;


% Find First and Second Law Efficiencies (eta_I, eta_II)
% ------------------------------------------------------
eta_I  = Wdot ./ (-dHdot);
eta_II = Wdot ./ (-dGdot);


specs.mol_air =         mol_air;
specs.mol_o2_react =    mol_o2_react;
specs.mol_n2 =          mol_n2;
% TODO update Spec to accomodate inlet water? (PLEASE CHECK KENDALL)

end