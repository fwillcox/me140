% findEtas.m
% Calculates the stoichiometry for a PEM fuel cell, dry air and dry h2
% and returns the first and second law efficiencies

% SOURCE: Lecture 8, Slides 21 & 22
% ASSUME:
% (i) 1 mol of H2                    
% (ii) isothermal & isobaric      % TODO: CHECK THESE ASSUMPTIONS!

function [eta_I, eta_II, Idot] = findEtas(T,Wdot)


global PERMIN_TO_PERSEC PERHR_TO_PERSEC G_PER_KG LHV F N_TO_O SCF_TO_MOLS C_TO_K
defineGlobals();
% Constants                                                  % Engineering Air Molar Mass Ratio of Nitrogen to Oxygen
specs = Spec();                                                     % class initialization

% Molar Masses
MM_h = 1.00794;
MM_o = 15.9994;
MM_n = 14.0067;
MM_h2o = 2*MM_h + MM_o;
MM_air = 28.97;

% Pressures
Pfuel =  [2.9 2.9 3.1 3.3 3.30 3.20 3.00 3.0 3.1];                  % psi (gauge)
Ptotal =   [0.2 0.3 0.6 0.7 1.15 1.25 1.35 1.3 1.5];                % psi (gauge), pressure of combined air and H2O after humidifier

% Mass Flow Rates (TODO: check what units the mdots should be in)
mdot_h2_scf =  [2.50 6.20 10.5 14.3 18.2 22.0 24.6 25.0 26.1];      % scf/hr (standard cubic feet/hour)
mdot_h2 = mdot_h2_scf .* SCF_TO_MOLS * PERHR_TO_PERSEC;             % mols/s
mdot_air_scf = [0.75 1.10 1.45 1.81 2.55 3.10 3.30 3.25 3.40];      % scf/min
mdot_air = mdot_air_scf * SCF_TO_MOLS * PERMIN_TO_PERSEC;           % mols/s
mdot_h2o = 40;                                                      % g/s (ASSUMPTION for Part 3) % TODO: convert this

% Excess Air Coefficient (lambda)
mol_h2 = 1;
mol_h2o_prod = mol_h2;
AF = mdot_air./mdot_h2;
mair = 0.5*(1+N_TO_O)*MM_air;
mfuel = mol_h2*(2*MM_h);
AFs = mair/mfuel;
lambda = AF./AFs;
mol_air = mair ./ (MM_air / G_PER_KG);

% Fraction of Liquid H2O (alpha)
% ASSUME: relative humidity = 100% = 1 --> alpha = Psat/Pair;
Psat = PsatW(T);
alpha = Ptotal./Psat; 

% Mols of Each Species
% --------------------
% Both (N2 stays constant)
mol_n2 = mol_air*N_TO_O/(1+N_TO_O);

% Reactants
mol_h2 = 1;                                                         % ASSUME: 1 mol of H2
mol_air = (1+N_TO_O)*lambda/2*mol_h2;
mol_o2_react = mol_air/(1+N_TO_O);
mol_total_react = mol_o2_react + mol_n2 + alpha;

% Products
mol_o2_prod = 0.5*(lambda-1).*mol_h2;                               % TODO: double check o2prod. Should it be *mol_h2?
mol_h2o_prod = mol_h2 + alpha;
beta = mol_h2o_prod - alpha;


% Mol Fractions of Each Species
% ------------------------------
y_o2_react = mol_o2_react ./mol_total_react;
y_n2_react = mol_n2       ./mol_total_react;
y_h2o_react = alpha       ./mol_total_react;

% Vapor Fraction of H2O 
% ---------------------
Pv_guess = Ptotal.*(beta./(beta + 0.5.*(lambda-1) +0.5.*lambda.*N_TO_O ));
if Pv_guess < Psat
    % All H2O is vapor (beta = 1)
    mol_h2ovap = beta;
    mol_h2oliq = 0;
else 
    % Some H2O is vapor, some liquid (beta not = 1)
    Pv_h2o = Psat;                                                  % LET: Pv = Psat to solve for beta
    y_h2o_prod = Pv_h2o./Ptotal; 
    beta = ((1+N_TO_O) .* y_h2o_prod)./ (1 - y_h2o_prod);
    mol_h2ovap = beta; 
    mol_h2oliq = mol_h2o_prod - mol_h2ovap;
end
pctVap = mol_h2ovap./(mol_h2o_prod);                                     % Percent of H2O that is vapor (by mass)


% Mol Fractions of Products
% --------------------------
mol_total_prod  = mol_o2_prod + mol_n2 + mol_h2ovap;
y_h2ovap  = mol_h2ovap  ./ mol_total_prod;
y_o2_prod = mol_o2_prod ./ mol_total_prod;
y_n2_prod = mol_n2      ./ mol_total_prod;


% Change in Gibbs Free Energy (dG)
% --------------------------------
% Reactants
Gdot_react = gEng(T, Pfuel,             'h2',mol_h2)      .* mdot_h2 ...
           + gEng(T, Ptotal .* y_o2_react,'o2',mol_o2_react).* mdot_air ...
           + gEng(T, Ptotal .* y_n2_react,'n2',mol_n2)      .* mdot_air;
if(alpha ~= 0 ) 
    Gdot_react = Gdot_react + gEng(T,Ptotal.*y_h2o_react,'h2ovap',alpha).* mdot_h2o; 
end

% Products
Gdot_prod = gEng(T, Ptotal.*y_h2ovap,   'h2ovap', mol_h2ovap) .* mdot_h2o... %% THIS WAS THE PART 2 BUG FROM PROJECT #4!
          + gEng(T, Ptotal,            'h2o',    mol_h2oliq) .* mdot_h2o...
          + gEng(T, Ptotal.*y_o2_prod, 'o2',     mol_o2_prod).* mdot_air...   
          + gEng(T, Ptotal.*y_n2_prod, 'n2',     mol_n2)     .* mdot_air;

dGdot = Gdot_prod - Gdot_react; 


% Change in Enthalpy (dH)
% -----------------------
% Reactants
Hdot_react = hEng(T,'h2', mol_h2)       .* mdot_h2... 
           + hEng(T, 'o2', mol_o2_react).* mdot_air...
           + hEng(T, 'n2', mol_n2)      .* mdot_air;
if(alpha ~= 0 )
    Hdot_react = Hdot_react + hEng(T,'h2ovap',alpha) .* mdot_h2o; 
end

% Products
Hdot_prod = hEng(T,'h2ovap', mol_h2ovap) .* mdot_h2o...
          + hEng(T,'h2o',    mol_h2oliq) .* mdot_h2o...
          + hEng(T,'o2',     mol_o2_prod).* mdot_air...
          + hEng(T,'n2',     mol_n2)     .* mdot_air;
      
dHdot = Hdot_prod - Hdot_react;


% Rate of Work (Wdot)
% -------------------
Ndot = mol_h2 .*mdot_h2;         
Idot = -Wdot - dGdot;              


% Find First and Second Law Efficiencies (eta_I, eta_II)
% ------------------------------------------------------
eta_I  = Wdot ./ (-dHdot);
eta_II = Wdot ./ (-dGdot);


specs.mol_air =      mol_air;
specs.mol_o2_react = mol_o2_react;
specs.mol_n2 =       mol_n2;
% TODO update Spec to accomodate inlet water? (PLEASE CHECK KENDALL)

end