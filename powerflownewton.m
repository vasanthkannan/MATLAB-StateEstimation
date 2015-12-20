%function [systemMWloss] = powerflow
print = 3;
% Power Flow 

% Def of Data 
% iter:			Iteration number
% Maxiter:		Maximum number of iteration
% powerflow_tolerance:   Maximum error allowed
% numbus:		Number of nodes in the system
% numline:		Number of lines in the system
% frombus(i), tobus(i):	Nodes from and to bus numbers for line i 
% R(i):			Resistance line i
% X(i):			Reactance line i
% Bcap(i)      Total line charging for line i
% Bustype(i):	Kind of bus i (l = PQ, g = PV, s = refbus)
% v(i):			Voltage of the node i
% Psched(i):	Active power inyected in the node i
% Qsched(i):	Reactive power inyected in the node i
% Vsched(i):	Voltage schedule for the node i

% Power Flow starts here

% Init voltage and angle Variables
v=ones(1,numbus);
angle=zeros(numbus,1);

%Build the Psched and Qsched tables
totload = 0.0;
for i = 1:numbus
   Psched(i)= -Pload(i);
   Qsched(i)= -Qload(i);
   totload = totload + Pload(i);
end

for g = 1:numgen
   i = genbus(g);
   Psched(i) = Psched(i) + Pgen(g);
   Qsched(i) = Qsched(i) + Qgen(g);
   v(i)=Vsched(g);
end

if print >= 3
   fprintf(' %s\n','Iter   MAXDP  MAXDPbus   MAXDQ  MAXDQbus')
end

Iter = 0;
converge = 0;


while 1
   
Iter = Iter + 1;   
MAXDP = 0.0;
MAXDPbus = 0;
MAXDQ = 0.0;
MAXDQbus = 0;
 
% Calculation of Pinj and Qinj using the v and angles values.
for i = 1: numbus
  if Bustype(i) == 'L' | Bustype(i) == 'G'
     Pinj(i)=0;
     Qinj(i)=0;
     for j = 1:numbus
        if i == j
         Pinj(i) = Pinj(i) + G(i,i)*v(i)^2;
         Qinj(i) = Qinj(i) - B(i,i)*v(i)^2;
        else         
         Pinj(i) = Pinj(i) + v(i)*v(j)*(G(i,j)*cos(angle(i)-angle(j))+B(i,j)*sin(angle(i)-angle(j))); 
         Qinj(i) = Qinj(i) + v(i)*v(j)*(G(i,j)*sin(angle(i)-angle(j))-B(i,j)*cos(angle(i)-angle(j)));
        end
     end 
     DelP(i)=Psched(i)-Pinj(i);
     DelQ(i) = 0.0;
     if Bustype(i) == 'L'
        DelQ(i)=Qsched(i)-Qinj(i);
     end
     if abs(DelP(i)) > MAXDP
      MAXDP = abs(DelP(i));
      MAXDPbus = i;
     end
     if abs(DelQ(i)) > MAXDQ
      MAXDQ = abs(DelQ(i));
      MAXDQbus = i;
     end
  elseif Bustype(i) == 'S'
     Pinj(i) = 0;
     Qinj(i) = 0;
     DelP(i) = 0;
     DelQ(i) = 0;
  end
end

% Store ones on the diagonal of the Jacobian

numbus2 = 2 * numbus;
for i = 1:numbus2
   for j = 1:numbus2
      J(i,j) = 0.0;
      if i == j 
         J(i,j) = 1.0;
      end
   end
end

% Build the Jacobian
% Note, the Jacobian is built with P rows(prow) and Q rows (qrow), and with Theta Comumns (thcol) 
% and V rows (vrow). 
% The Jacobian is built one Prow and One Q row at a time. 

for i = 1: numbus
   prow = 2*i -1;
   qrow = 2*i;
 if Bustype(i) == 'L' | Bustype(i) == 'G'   
%   Build JPTh terms
    for j = 1: numbus
    thcol = 2*j -1;
    vcol =  2*j;
      if Bustype(j) == 'L' | Bustype(j) == 'G'
        if i == j    
          J(prow,thcol) = -Qinj(i) - B(i,i)*v(i)^2; 
        else
          J(prow,thcol) = v(i)*v(j)*(G(i,j)*sin(angle(i)-angle(j)) - B(i,j)*cos(angle(i)-angle(j)));
        end
      end
    end % end build JPth terms loop
%   Build matrix JPV
    for j = 1: numbus
    thcol = 2*j -1;
    vcol =  2*j;
      if Bustype(j) == 'L'
        if i == j    
          J(prow,vcol) = Pinj(i) + G(i,i)*v(i)^2; 
        else
          J(prow,vcol) = v(i)*v(j)*(G(i,j)*cos(angle(i)-angle(j))+B(i,j)*sin(angle(i)-angle(j)));
        end
      end
    end % end build JPV term loop
%   Build JQTh terms
   if Bustype(i) == 'L'
    for j = 1: numbus
    thcol = 2*j -1;
    vcol =  2*j;
     if Bustype(j) == 'L' | Bustype(j) == 'G'
        if i == j    
          J(qrow,thcol) = Pinj(i) - G(i,i)*v(i)^2;
        else
          J(qrow,thcol)= -v(i)*v(j)*(G(i,j)*cos(angle(i)-angle(j))+B(i,j)*sin(angle(i)-angle(j)));
        end
      end
     end %  end Build JQTh terms loop
    end  % end if   
%   Build JQV  terms
   if Bustype(i) == 'L'
    for j = 1: numbus
    thcol = 2*j -1;
    vcol =  2*j;
      if Bustype(j) == 'L'
        if i == j    
          J(qrow,vcol) = Qinj(i) - B(i,i)*v(i)^2;
        else
          J(qrow,vcol) = v(i)*v(j)*(G(i,j)*sin(angle(i)-angle(j))-B(i,j)*cos(angle(i)-angle(j)));
        end
      end
     end% end build JQV terms loop
   end % end if statement
   
 end % end if statement on L or G
   DelPQ(prow) = DelP(i);
   DelPQ(qrow) = DelQ(i);  

end % end master bus loop building Jacobian
if converge
    break    
end

% Calculation of the delta voltages and delta angles
%DelTHV=inv(J)*DelPQ';
DelTHV = J\DelPQ';

% Assign difference calculated in last step
for i = 1: numbus
   if Bustype(i) == 'L' | Bustype(i) == 'G'
     throw = 2*i-1;
     angle(i) = angle(i) + DelTHV(throw);
     if Bustype(i) == 'L'
        vrow = 2*i;
        v(i) = v(i)*(1 + DelTHV(vrow));
     end
   end
end;

% Print and save result from last iteration
if print >= 3
   fprintf(' %2d %10.6f %6d %10.6f %5d\n',Iter, MAXDP, MAXDPbus, MAXDQ, MAXDQbus)
end

if MAXDP < powerflow_tolerance
   if MAXDQ < powerflow_tolerance
      converge = 1;
      for g = 1:numgen
         Qgen(g) = Qinj(genbus(g)) + Qload(genbus(g));
         if Qgen(g) > Qmax(g)
          Bustype(genbus(g)) = 'L';
          Qsched(genbus(g)) = Qmax(g) - Qload(genbus(g));
          converge = 0;
         end
      end
   end
end
   
if Iter > Maxiter
  'Maximun number of iteration reached - The Power Flow has failed'
  converge = 1;
end 


end % end while loop

%V=v'; 
%V(find(Bustype=='G'))

% Display results 
% Calculate the net bus injections 
for i = 1: numbus
  Pinj(i)=0;
  Qinj(i)=0;
  for j = 1:numbus
     if i == j
      Pinj(i) = Pinj(i) + G(i,i)*v(i)^2;
      Qinj(i) = Qinj(i) - B(i,i)*v(i)^2;
     else         
      Pinj(i) = Pinj(i) + v(i)*v(j)*(G(i,j)*cos(angle(i)-angle(j))+B(i,j)*sin(angle(i)-angle(j))); 
      Qinj(i) = Qinj(i) + v(i)*v(j)*(G(i,j)*sin(angle(i)-angle(j))-B(i,j)*cos(angle(i)-angle(j)));
     end
  end
end 

totPload = 0.0;
totQload = 0.0;
totPinj = 0.;
totQinj = 0.;
for i = 1:numbus
   totPload = totPload + Pload(i)* baseMVA;
   totQload = totQload + Qload(i)* baseMVA;
   totPinj = totPinj + Pinj(i)* baseMVA;
   totQinj = totQinj + Qinj(i)* baseMVA;
end
totPload;
totQload;
totPinj;
totQinj;

totPgen = 0.0;
totQgen = 0.0;
for g = 1:numgen
   Pgen(g) = (Pinj(genbus(g)) + Pload(genbus(g)));
   Qgen(g) = (Qinj(genbus(g)) + Qload(genbus(g)));
   totPgen = totPgen + Pgen(g)* baseMVA;
   totQgen = totQgen + Qgen(g)* baseMVA;
end
totMWgen = totPgen;
totMWload = totPload ;
totMVARload = totQload ;

% Place voltage mag and angle into a complex vector
vr = v.*cos(angle') + sqrt(-1)*v.*sin(angle');

systemMWloss = 0.0;
for j = 1:numline
   Iline = (vr(frombus(j)) - vr(tobus(j)))/(R(j)+sqrt(-1)*X(j));
   lineMWloss = ( abs(Iline)^2 * R(j) )* baseMVA;
   systemMWloss = systemMWloss + lineMWloss;
end

Jinv = inv(J);

if print >= 1
fprintf(' %s \n',' ');
fprintf(' %s  %10.3f %s%10.3f \n','Power Flow with Total Pgen = ',totMWgen, ' Total Qgen = ',totQgen);
fprintf(' %s  %10.3f %s%10.3f\n', '                Total PLoad= ',totMWload,' Total Qload= ',totMVARload);
fprintf(' %s  %10.3f \n','           Total MW Losses = ',systemMWloss);
fprintf(' %s \n',' ');
end % end of print option print>1

for i = 1:numbus
   busgen(i) = 0;
end

for g = 1:numgen
      busgen(genbus(g))=g;
      PgenMW(g) = Pgen(g)*baseMVA;
      QgenMW(g) = Qgen(g)*baseMVA;
end

for i = 1:numbus
   PloadMW(i) = Pload(i)*baseMVA;
   QloadMW(i) = Qload(i)*baseMVA;

end

for i = 1:numbus
   Vmag(i) = v(i);
   angle_deg(i) = angle(i)*180/pi;
end

if print >= 2
fprintf(' %s \n',' ');
fprintf(' %s \n','Bus  Vmag     angle     Pgen     Qgen     Pload   Qload  To Bus    Pline   Qline ');
for i = 1:numbus
   if i >9
      space = '';
   else
      space = ' ';
   end
   
   if busgen(i) > 0
      fprintf(' %s%d  %5.3f   %6.3f %8.3f %8.3f %8.3f %8.3f \n',space,i,Vmag(i),angle_deg(i),PgenMW(busgen(i)),QgenMW(busgen(i)),PloadMW(i),QloadMW(i));
   else
      fprintf(' %s%d  %5.3f   %6.3f  %s%8.3f %8.3f \n',space,i,Vmag(i), angle_deg(i),'                 ',PloadMW(i),QloadMW(i)); 
   end
   
   for j = 1:numline
      if frombus(j) == i;
          if tobus(j) >9
             space = '';
          else
             space = ' ';
          end
          Iline = (vr(frombus(j)) - vr(tobus(j)))/(R(j)+sqrt(-1)*X(j));
          Iline = Iline + vr(frombus(j))*(0.0 + sqrt(-1)*Bcap(j)*0.5);
          Sline = vr(frombus(j)) * conj(Iline);
          Pline = real(Sline)*baseMVA;
          Qline = imag(Sline)*baseMVA;
          fprintf(' %s%s  %d  %8.3f %8.3f  \n',space,'                                                          ',tobus(j),Pline, Qline);
       elseif tobus(j) == i
          if frombus(j) >9
             space = '';
          else
             space = ' ';
          end
          Iline = (vr(tobus(j)) - vr(frombus(j)))/(R(j)+sqrt(-1)*X(j));
          Iline = Iline + vr(tobus(j))*(0.0 + sqrt(-1)*Bcap(j)*0.5);
          Sline = vr(tobus(j)) * conj(Iline);
          Pline = real(Sline)*baseMVA;
          Qline = imag(Sline)*baseMVA;
          fprintf(' %s%s  %d  %8.3f %8.3f  \n',space,'                                                          ',frombus(j),Pline, Qline);
       end
    end
 end
 
end % end of print option print>1

PowerFlowInputData.baseMVA      = baseMVA;
PowerFlowInputData.Maxiter      = Maxiter;
PowerFlowInputData.powerflow_tolerance    = powerflow_tolerance;
PowerFlowInputData.numbus       = numbus;
PowerFlowInputData.numline      = numline;
PowerFlowInputData.numgen       = numgen;
PowerFlowInputData.numarea      = numarea;
PowerFlowInputData.refbus       = refbus;
PowerFlowInputData.frombus      = frombus;
PowerFlowInputData.tobus        = tobus;
PowerFlowInputData.R            = R;
PowerFlowInputData.X            = X;
PowerFlowInputData.Bcap         = Bcap;
PowerFlowInputData.Bustype      = Bustype;
PowerFlowInputData.Psched       = Psched;
PowerFlowInputData.Qsched       = Qsched;
PowerFlowInputData.Vsched       = Vsched;
PowerFlowInputData.Y            = Y;
PowerFlowInputData.G            = G;
PowerFlowInputData.B            = B;

PowerFlowSolution.Pinj	   = Pinj;
PowerFlowSolution.Qinj	   = Qinj;
PowerFlowSolution.vr       = vr;
PowerFlowSolution.Vmag     = Vmag;
PowerFlowSolution.Theta    = angle;
PowerFlowSolution.Bustype  = Bustype;

save('PowerFlowOutput', 'PowerFlowSolution','PowerFlowInputData');


