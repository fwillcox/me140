%gEng.m , wrapper for energyF
% 4-22-16 created by Jon Renslo
function G = gEng(T,P,spec,mol)
%         for i = 1:length(T)
%             if length(P) == 1
%                 j = 1;
%             else
%                 j = i;
%             end
%             if length(mol) == 1
%                 k = 1;
%             else
%                 k = i;
%             end
% 
%             if(nargin == 3)     % Note: nargin returns number of function input arguments
%                 mol = 1;
%             end
%             a = energyF(T(i),P(j),spec,mol(k));
%             G(i) = a.G;
%         end
    a = energyF(T,P,spec,mol);
    G = a.G;
end