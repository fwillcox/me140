% findEtas.m
% Calculates the stoichiometry for a PEM fuel cell, dry air and dry h2
% and returns the first and second law efficiencies

% SOURCE: Lecture 8, Slides 21 & 22
% ASSUME:
% (i) 1 mol of H2
% (ii) isothermal & isobaric      % TODO: CHECK THESE ASSUMPTIONS!

function [eta_I, eta_II, Idot,lambda] = findEtas(mdot_wetair_react, mdot_fuel, Ptotal, Pfuel, T, Wdot)

global PERMIN_TO_PERSEC PERHR_TO_PERSEC G_PER_KG LHV F N_TO_O SCF_TO_MOLS ...
    C_TO_K PSI_TO_PA MM_h MM_h2 MM_o MM_n MM_h2o MM_air PATM
defineGlobals();
% Constants                                                  % Engineering Air Molar Mass Ratio of Nitrogen to Oxygen
specs = Spec();                                                     % class initialization

% Find Excess Air Coefficient (lambda) by assuming 1mol h2
mol_h2 = 1; % ASSUME: 1 mol of H2
mairstoich = 0.5*(1+N_TO_O)*MM_air;
mfuelstoich = mol_h2*(MM_h2);
AFs = mairstoich/mfuelstoich;

AF = mdot_wetair_react./mdot_fuel;
lambda = AF./AFs;

% Mols of Each Species
% --------------------
% Both (N2 stays constant)
moldot_h2 = mdot_fuel / (MM_h2 / G_PER_KG);
moldot_air_dry_react = (1+N_TO_O) .* lambda / 2 .* moldot_h2;
moldot_n2 = moldot_air_dry_react .* N_TO_O/(1+N_TO_O);

% Mols of H2O in Reactants (alpha)
% ASSUME: relative humidity = 100% = 1 --> alpha = Psat/Pair;
Psat = PsatW(T);
alpha = Ptotal./Psat .* moldot_h2;

% Reactants
moldot_o2_react = moldot_air_dry_react./(1+N_TO_O);
moldot_total_react = moldot_o2_react + moldot_n2 + alpha;

% Products
moldot_o2_prod = 0.5*(lambda-1).*moldot_h2;
moldot_h2o_prod = moldot_h2 + alpha;
beta = moldot_h2o_prod;                 %assume its all vapor initially

% Mol Fractions of Each Species
% ------------------------------
y_o2_react = moldot_o2_react ./moldot_total_react;
y_n2_react = moldot_n2       ./moldot_total_react;
y_h2o_react = alpha       ./moldot_total_react;

% Vapor Fraction of H2O (correcting beta and gamma)
% ---------------------
Pv_guess = Ptotal.*(beta./(beta + 0.5.*(lambda-1) +0.5.*lambda.*N_TO_O ));
if Pv_guess < Psat
    % All H2O is vapor (beta = 1)
    moldot_h2ovap = beta;
    moldot_h2oliq = 0;
else
    % Some H2O is vapor, some liquid (beta not = 1)
    Pv_h2o = Psat;                                                  % LET: Pv = Psat to solve for beta
    y_h2o_prod = Pv_h2o./Ptotal;
    beta = ((1+N_TO_O) .* y_h2o_prod)./ (1 - y_h2o_prod);
    moldot_h2ovap = beta;
    moldot_h2oliq = moldot_h2o_prod - moldot_h2ovap;
end
pctVap = moldot_h2ovap./(moldot_h2o_prod);                                     % Percent of H2O that is vapor (by mass)


% Mol Fractions of Products
% --------------------------
moldot_total_prod  = moldot_o2_prod + moldot_n2 + moldot_h2ovap;
y_h2ovap_prod  = moldot_h2ovap  ./ moldot_total_prod;
y_o2_prod = moldot_o2_prod ./ moldot_total_prod;
y_n2_prod = moldot_n2      ./ moldot_total_prod;


% Change in Gibbs Free Energy (dG)
% --------------------------------
% Reactants
Gdot_react = gEng(T, Pfuel,             'h2'  ,moldot_h2)       ...   % [J*s^-1]
    + gEng(T, Ptotal .* y_o2_react,'o2',moldot_o2_react) ...
    + gEng(T, Ptotal .* y_n2_react,'n2',moldot_n2)      ;
if(alpha ~= 0 )
    Gdot_react = Gdot_react + gEng(T,Ptotal.*y_h2o_react,'h2ovap',alpha);
end

% Products
Gdot_prod = gEng(T, Ptotal.*y_h2ovap_prod,   'h2ovap', moldot_h2ovap) ...
    + gEng(T, Ptotal,            'h2o',    moldot_h2oliq)  ...
    + gEng(T, Ptotal.*y_o2_prod, 'o2',     moldot_o2_prod)...
    + gEng(T, Ptotal.*y_n2_prod, 'n2',     moldot_n2)   ;

dGdot = Gdot_prod - Gdot_react;

%check mol fractions
y_o2_react + y_n2_react + y_h2o_react;
y_h2ovap_prod + y_o2_prod + y_n2_prod;
%check mass balance
masscreation = moldot_h2*MM_h2 + moldot_o2_react*MM_o*2 + moldot_n2*MM_n*2 ...
    + alpha*MM_h2o - ...
    ( (moldot_h2ovap + moldot_h2oliq) * MM_h2o + moldot_o2_prod * MM_o*2 ...
    + moldot_n2*MM_n*2);


% Change in Enthalpy (dH)
% -----------------------
% Reactants
Hdot_react = hEng(T,'h2', moldot_h2)       ... % [J*s^-1]
    + hEng(T, 'o2', moldot_o2_react)...
    + hEng(T, 'n2', moldot_n2)      ;
if(alpha ~= 0 )
    Hdot_react = Hdot_react + hEng(T,'h2ovap',alpha);
end

% Products
Hdot_prod = hEng(T,'h2ovap', moldot_h2ovap) ... % [J*s^-1]
    + hEng(T,'h2o',    moldot_h2oliq) ...
    + hEng(T,'o2',     moldot_o2_prod)...
    + hEng(T,'n2',     moldot_n2)     ;

dHdot = Hdot_prod - Hdot_react;


% Rate of Work (Wdot)
% -------------------
Idot = -Wdot - dGdot;

% Find First and Second Law Efficiencies (eta_I, eta_II)
% ------------------------------------------------------
eta_I  = Wdot ./ (LHV*mdot_fuel);       % Watts / (J/kg * kg/s)
eta_II = (Wdot) ./ (-dGdot);  % Watts / Watts


specs.mol_air =      moldot_air_dry_react;
specs.mol_o2_react = moldot_o2_react;
specs.mol_n2 =       moldot_n2;

end