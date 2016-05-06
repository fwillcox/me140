function  defineGlobals( )

global PERMIN_TO_PERSEC PERHR_TO_PERSEC G_PER_KG LHV F N_TO_O SCF_TO_MOLS ...
    C_TO_K PSI_TO_PA MM_h MM_h2 MM_o MM_n MM_h2o MM_ch4 MM_air PATM HORSEPOWER_TO_W  MM_c
% Temperature
C_TO_K = 273.15;

% Mass
G_PER_KG = 1000;        

% Time
PERHR_TO_PERSEC = 1/(60*60);
PERMIN_TO_PERSEC = 1/60;

% Volume
SCF_TO_MOLS = 1.19804;  % standard cubic feet to moles

% Pressure
PSI_TO_PA = 6894.76;

% Power
HORSEPOWER_TO_W = 745.7;        % horsepower to watts

% Useful Constants
LHV = 120.0*10^6;       % J/kg,  Lower Heating Value H2  
F = 96485;              % C/(mol of e-), Faraday's Constant
N_TO_O = 79/21;  

% Molar Masses
MM_h = 1.00794; % g/mol
MM_h2 = 2* MM_h;
MM_o = 15.9994;
MM_n = 14.0067;
MM_c = 12.01;
MM_h2o = 2*MM_h + MM_o;
MM_air = 28.97;
MM_ch4 = MM_c + 4*MM_h;

PATM = 101.3e3; %pascals


end

