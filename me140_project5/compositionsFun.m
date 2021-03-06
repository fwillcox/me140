% Returns an array with the composition of each species at a given kp.  
% Species index order: CO, H20, CO2, H2
function comps = compositionsFun(kp)
     syms nco nh2o nco2 nh2
     eqs = [  1  == nco2   + nco;...          carbon atom balance  %POTENTIAL ERROR: shouldn't this be 2, not 1? NATASHA: THINK 1 WORKS BUT SOMEONE ELSE SHOULD ALSO CHECK
         3  == nco2*2 + nco + nh2o; ...  oxygen atom balance
         10  == nh2*2   + nh2o*2;...      hydrogen atom balance
         (nco2.*nh2)./(nco.*nh2o) == kp];  %Nernst atom balance % Tin(1) is the first temp in Tin vector, which is Tin(reformer)  
     assume([nco,nh2o,nco2,nh2],'real'); assumeAlso([nco,nh2o,nco2,nh2] > 0); assumeAlso([nco,nh2o,nco2,nh2] < 20)
     [a,b,c,d] = vpasolve(eqs,[nco,nh2o,nco2,nh2],[1,1,1,1]);
     sol = [a,b,c,d];
     sizeSent = isequal([4 1],size(sol'));
     comps(:,1) = double(sol)';
end
