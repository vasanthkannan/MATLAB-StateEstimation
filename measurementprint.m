if print >= 2
fprintf(' %s \n',' ');
fprintf(' %s \n','               Base Case Value    Measured Value ');
fprintf(' %s \n','              _________________  _________________');
fprintf(' %s \n','Measurement   kV    MW    MVAR   kV    MW    MVAR ');
for i = 1:numbus
   if i >9
      space = '';
   else
      space = ' ';
   end
   
   if numVmeas > 0
       for j = 1:numVmeas
           if Vmeasbus(j) == i 
            fprintf(' %s  %5.1f   %5.1 \n',voltage_meas_NAMES(j,:),Vmeasvalueinput(j),Vmeasvalue(j));
           end
       end
   end
   
%    for j = 1:numline
%       if frombus(j) == i;
%           if tobus(j) >9
%              space = '';
%           else
%              space = ' ';
%           end
%           Iline = (vr(frombus(j)) - vr(tobus(j)))/(R(j)+sqrt(-1)*X(j));
%           Iline = Iline + vr(frombus(j))*(0.0 + sqrt(-1)*Bcap(j)*0.5);
%           Sline = vr(frombus(j)) * conj(Iline);
%           Pline = real(Sline)*baseMVA;
%           Qline = imag(Sline)*baseMVA;
%           fprintf(' %s%s  %d  %8.3f %8.3f  \n',space,'                                                          ',tobus(j),Pline, Qline);
%        elseif tobus(j) == i
%           if frombus(j) >9
%              space = '';
%           else
%              space = ' ';
%           end
%           Iline = (vr(tobus(j)) - vr(frombus(j)))/(R(j)+sqrt(-1)*X(j));
%           Iline = Iline + vr(tobus(j))*(0.0 + sqrt(-1)*Bcap(j)*0.5);
%           Sline = vr(tobus(j)) * conj(Iline);
%           Pline = real(Sline)*baseMVA;
%           Qline = imag(Sline)*baseMVA;
%           fprintf(' %s%s  %d  %8.3f %8.3f  \n',space,'                                                          ',frombus(j),Pline, Qline);
%        end
%     end
%  end
 
end % end of print option print>1
