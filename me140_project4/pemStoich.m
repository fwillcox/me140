%pemStoich.m
%4-22-16 Created Jon Renslo
function [ mixVec ] = pemStoich( lambda, m_h2)
%PEMSTOICH 
% does stoichiometry calculation for a PEM H2 and air fuel cell
% mH2 in kg
% mixVec = [ water vapor, liquid water, o2, n2] in moles of product formed
MM_h2 = 2.02; %g/mol
G_TO_KG = 1000;

N_h2 = m_h2/MM_h2*G_TO_KG;

%TODO based on exit conditions, must calculate liquid/gas water phase for products

mixVec = [0];


end