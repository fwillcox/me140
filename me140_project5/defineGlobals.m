function  defineGlobals( )

global PERMIN_TO_PERSEC PERHR_TO_PERSEC G_PER_KG LHV F N_TO_O SCF_TO_MOLS C_TO_K PSI_TO_PA
% Temperature
C_TO_K = 273.15;

% Mass
G_PER_KG = 1000;

% Time
PERHR_TO_PERSEC = 1/(60*60);
PERMIN_TO_PERSEC = 1/60;

% Volume
SCF_TO_MOLS = 1.19804;

% Pressure
PSI_TO_PA = 6894.76;

% Useful Constants
LHV = 120.0*10^6;     % J/kg,  Lower Heating Value H2  
F = 96485;            % C/(mol of e-), Faraday's Constant
N_TO_O = 79/21;  



end

