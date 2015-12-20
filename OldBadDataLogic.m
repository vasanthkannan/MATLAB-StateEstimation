%    fprintf(' %s \n','');
%    fprintf(' %s \n','');
%    fprintf(' %s \n','');
%    fprintf(' %s \n','-------------------');
%    fprintf(' %s \n','Bad Data Analysis');
%    fprintf(' %s \n','-----------------------------------------------------------------------');
%    fprintf(' %s \n','Measurement  Measurement    Caculated  Difference   Covariance   y_norm ');
%    fprintf(' %s \n','Number       Value          Value                   of error            ');
%    fprintf(' %s \n','-----------------------------------------------------------------------');
%    fprintf(' %s\n','Pflow measurements:');
%    for i=1:numFmeas
%        fprintf('%5d       %8.3f       %8.3f     %8.3f     %8.3f    %8.3f\n', i, FmeasPvalue(i)*baseMVA,  FcalcPvalue(i)*baseMVA, (FmeasPvalue(i)-FcalcPvalue(i))*baseMVA,sigmaPflow(i),y1_norm(i)*baseMVA);
%    end
%    fprintf(' %s\n','Qflow measurements:');
%    for i=1:numFmeas
%        fprintf('%5d       %8.3f       %8.3f     %8.3f     %8.3f    %8.3f\n', i, FmeasQvalue(i)*baseMVA,  FcalcQvalue(i)*baseMVA, (FmeasQvalue(i)-FcalcQvalue(i))*baseMVA,sigmaQflow(i),y2_norm(i)*baseMVA);
%    end 
%    fprintf(' %s\n','Pinjection measurements:');
%    for i=1:numImeas
%        fprintf('%5d       %8.3f       %8.3f     %8.3f     %8.3f    %8.3f\n', Imeasbus(i), ImeasPvalue(i)*baseMVA,  IcalcPvalue(i)*baseMVA, (ImeasPvalue(i)-IcalcPvalue(i))*baseMVA,sigmaPinjec(i),y3_norm(i)*baseMVA);
%    end
%    fprintf(' %s\n','Qinjection measurements:');
%    for i=1:numImeas
%        fprintf('%5d       %8.3f       %8.3f     %8.3f     %8.3f    %8.3f\n', Imeasbus(i), ImeasQvalue(i)*baseMVA,  IcalcQvalue(i)*baseMVA, (ImeasQvalue(i)-IcalcQvalue(i))*baseMVA,sigmaQinjec(i),y4_norm(i)*baseMVA);
%    end 
% %find the maximum y_norm(i), which is the prime suspect of error,both
% %the value and the index, putting them into Primesuspectvalue and
% %Primesuspectindex seperately.
% 
% 
% [Primesuspectvalue, Primesuspectindex]=max(abs(y_norm))
% if Primesuspectvalue<0.03
%     fprintf('There is no bad measurement anymore!\n');
%     break;
% end
% 
% % Remove the bad measurement and display the bad measurement index and
% % value.
% if Primesuspectindex<=numFmeas    % the removed measurements is a P flow measurement
%     fprintf('The bad measurement is a P flow measurement. %e\n');
%     fprintf('The bad measurement occurs at(the index of P flow measurements): %2d    %s', Primesuspectindex,'');
%     fprintf('Its value is:%e',Primesuspectvalue);
%     fprintf(' %s \n',' ');
%     numFmeas=numFmeas-1;    
%     if numFmeas~=0
%       for i=Primesuspectindex:(numFmeas)
%         FmeasPvalue(i)=FmeasPvalue(i+1);
%         FmeasPsigma(i)=FmeasPsigma(i+1);
%         Fmeasbranch(i)=Fmeasbranch(i+1);
%         Fmeasfrombus(i)=Fmeasfrombus(i+1);
%         Fmeastobus(i)=Fmeastobus(i+1);
%       end
%     end
% elseif Primesuspectindex>numFmeas&&Primesuspectindex<=(numFmeas+numFmeas)
%     fprintf('The bad measurement is a Q flow measurement. %e\n');
%     fprintf('The bad measurement occurs at (the index of Q flow measurements):  %2d    %s', (Primesuspectindex-numFmeas),'');
%     fprintf('Its value is:%e',Primesuspectvalue);
%     fprintf(' %s \n',' ');
%     numFmeas=numFmeas-1;
%     if numFmeas~=0
%       for i=(Primesuspectindex-numFmeas):(numFmeas)
%         FmeasQvalue(i)=FmeasQvalue(i+1);
%         FmeasQsigma(i)=FmeasQsigma(i+1);
%         Fmeasbranch(i)=Fmeasbranch(i+1);
%         Fmeasfrombus(i)=Fmeasfrombus(i+1);
%         Fmeastobus(i)=Fmeastobus(i+1);
%       end
%     end
% elseif Primesuspectindex>(numFmeas+numFmeas) && Primesuspectindex<=(numFmeas+numFmeas+numImeas)
%     badmeasurementbus=Primesuspectindex-numFmeas-numFmeas;
%     fprintf('The bad measurement is a P injection measurement. %e\n');
%     fprintf('The bad measurement occurs at bus: %2d    %s',  Imeasbus(badmeasurementbus),'');
%     fprintf('Its value is:%2d',Primesuspectvalue);
%     fprintf(' %s \n',' ');
%     numImeas=numImeas-1;
%     if numImeas~=0
%       for i=(Primesuspectindex-numFmeas-numFmeas):(numImeas)
%         ImeasPvalue(i)=ImeasPvalue(i+1);
%         ImeasPsigma(i)=ImeasPsigma(i+1);
%         Imeasbus(i)=Imeasbus(i+1);
%       end
%       Imeasbus=Imeasbus(1:numImeas);
%     end
% else
%     badmeasurementbus=Primesuspectindex-numFmeas-numFmeas-numImeas;
%     
% end
% end
%     fprintf('The bad measurement is a Q injection measurement. %e\n');
%     fprintf('The bad measurement occurs at bus: %2d    %s',  Imeasbus(badmeasurementbus), '');
%     fprintf('Its value is:%2d',Primesuspectvalue);
%     fprintf(' %s \n',' ');
%     numImeas=numImeas-1;
%     if numImeas~=0
%       for i=(Primesuspectindex-numFmeas-numFmeas-numImeas):(numImeas)
%         ImeasQvalue(i)=ImeasQvalue(i+1);
%         ImeasQsigma(i)=ImeasQsigma(i+1);
%         Imeasbus(i)=Imeasbus(i+1);
%       end
%       Imeasbus=Imeasbus(1:numImeas);
%     end
% end
%       % Calculate all f(i),y_norm(i) values for Nm measurements.
%   for i=1:numbus
%      vr(i)=z.vmag(i)*(cos(z.vang(i))+sqrt(-1)*sin(z.vang(i)));
%   end
%   for i = 1:numFmeas
%     a=Fmeasfrombus(i);
%     b=Fmeastobus(i);
% %         for j = 1:numline
% %           if frombus(j) == a & tobus(j)==b
% %               break
% %           elseif frombus(j) == b & tobus(j)==a
% %               break
% %           end
% %         end
%     Iline = (vr(a) - vr(b))/(R(Fmeasbranch(i))+sqrt(-1)*X(Fmeasbranch(i)));
%     Iline = Iline + vr(a)*(0.0 + sqrt(-1)*Bcap(Fmeasbranch(i))*0.5);
%     Sline = vr(a) * conj(Iline);
%     FcalcPvalue(i) = real(Sline);
%     y1_norm(i)=(FmeasPvalue(i)-FcalcPvalue(i))/sigmaPflow(i);
% 
%   end
%    y1_norm=y1_norm(1:numFmeas);
% 
% for i = 1:numFmeas
%     a=Fmeasfrombus(i);
%     b=Fmeastobus(i);
% %         for j = 1:numline
% %           if frombus(j) == a & tobus(j)==b
% %               break
% %           elseif frombus(j) == b & tobus(j)==a
% %               break
% %           end
% %         end
%     Iline = (vr(a) - vr(b))/(R(Fmeasbranch(i))+sqrt(-1)*X(Fmeasbranch(i)));
%     Iline = Iline + vr(a)*(0.0 + sqrt(-1)*Bcap(Fmeasbranch(i))*0.5);
%     Sline = vr(a) * conj(Iline);
%     FcalcQvalue(i)= imag(Sline);
%     %if FmeasQvalue(i)-FcalcQvalue(i)<0.00001
%      % y2_norm(i)=0;
%    %else
%       y2_norm(i)=(FmeasQvalue(i)-FcalcQvalue(i))/sigmaQflow(i);
% 
%    % end
% end
% y2_norm=y2_norm(1:numFmeas);
% 
% IcalcPvalue=zeros(numImeas,1);
% IcalcQvalue=zeros(numImeas,1);
% for i=1:numImeas
%     a = Imeasbus(i);
%     for b=1:numbus
%         Pinj =z.vmag(a)*z.vmag(b)*(cos(z.vang(a)-z.vang(b))*G(a,b)-sin(z.vang(a)-z.vang(b))*B(a,b));
%         IcalcPvalue(i)=IcalcPvalue(i)+Pinj;
%     end
%     y3_norm(i)=(ImeasPvalue(i)-IcalcPvalue(i))/sigmaPinjec(i);
% 
% end
%  y3_norm=y3_norm(1:numImeas);
% for i=1:numImeas
%     a = Imeasbus(i);
%     for b=1:numbus
%         Qinj = z.vmag(a)*z.vmag(b)*(sin(z.vang(a)-z.vang(b))*G(a,b)+cos(z.vang(a)-z.vang(b))*B(a,b));
%         IcalcQvalue(i)=IcalcQvalue(i)+Qinj;
%     end
%     y4_norm(i)=(ImeasQvalue(i)-IcalcQvalue(i))/sigmaQinjec(i);
% 
% end
% y4_norm=y4_norm(1:numImeas);
% % Assume that numFmeas, numFmeas never go to zero.
% if numImeas==0
%     y_norm=[y1_norm,y2_norm,y4_norm];
% elseif numImeas==0
%     y_norm=[y1_norm,y2_norm,y3_norm];
% elseif numImeas==0&&numImeas==0
%     y_norm=[y1_norm,y2_norm];
% else
%     y_norm=[y1_norm,y2_norm,y3_norm,y4_norm];
% end
% end       % end "Detect and Identify bad measurements" process
% 
% %end % end threshold iteration loop
% 
% %Recalculate y_norm
% for i = 1:numFmeas
% y1_norm(i)=(FmeasPvalue(i)-FcalcPvalue(i))/sigmaPflow(i);
% end
% for i=1:numFmeas
% y2_norm(i)=(FmeasQvalue(i)-FcalcQvalue(i))/sigmaQflow(i);
% end
% for i=1:numImeas
% y3_norm(i)=(ImeasPvalue(i)-IcalcPvalue(i))/sigmaPinjec(i);
% end
% for i=1:numImeas
% y4_norm(i)=(ImeasQvalue(i)-IcalcQvalue(i))/sigmaQinjec(i);
% end
% 
% fprintf(' %s \n',' ');
% fprintf ('The final value of J is %d\n' , J)
% fprintf ('The final value of TJ is %d\n' , TJ) 
%  % Display
%     fprintf(' %s \n',' ');
%     fprintf(' %s \n','-----------------------------------------------------------------------------------');
%     fprintf(' %s \n','From    To       Measured       Calculated        Difference   Covariance of       ');
%     fprintf(' %s \n','Bus     Bus      Value          Value                          Meansurement        ');
%     fprintf(' %s \n','-----------------------------------------------------------------------------------');
% 
%     fprintf(' %s \n','Pflow: ');
%     for i=1:numFmeas
%         fprintf(' %5d %5d    %8.3f        %8.3f        %8.3f       %8.3f  \n',Fmeasfrombus(i),Fmeastobus(i),FmeasPvalue(i)*baseMVA, FcalcPvalue(i)*baseMVA,(FmeasPvalue(i)-FcalcPvalue(i))*baseMVA,FmeasPsigma(i));
%     end
%     fprintf(' %s \n',' ');
%     fprintf(' %s \n',' ');
% 
%     fprintf(' %s \n','Qflow: ');
%     for i=1:numFmeas
%         fprintf(' %5d %5d    %8.3f        %8.3f        %8.3f       %8.3f  \n',Fmeasfrombus(i),Fmeastobus(i),FmeasQvalue(i)*baseMVA, FcalcQvalue(i)*baseMVA,(FmeasQvalue(i)-FcalcQvalue(i))*baseMVA,FmeasQsigma(i));
%     end
%     fprintf(' %s \n',' ');
%     fprintf(' %s \n',' ');
% 
% 
%     fprintf(' %s \n','Pinjection: ');
%     for i=1:numImeas
%         fprintf(' %5d           %8.3f        %8.3f       %8.3f       %8.3f     \n',Imeasbus(i),ImeasPvalue(i)*baseMVA, IcalcPvalue(i)*baseMVA,(ImeasPvalue(i)-IcalcPvalue(i))*baseMVA, ImeasPsigma(i));
%     end
% 
%     fprintf(' %s \n',' ');
%     fprintf(' %s \n',' ');
%     fprintf(' %s \n','Qinjection: ');
% 
%     for i=1:numImeas  
%         fprintf(' %5d           %8.3f        %8.3f       %8.3f       %8.3f     \n',Imeasbus(i),ImeasQvalue(i)*baseMVA, IcalcQvalue(i)*baseMVA,(ImeasQvalue(i)-IcalcQvalue(i))*baseMVA,ImeasQsigma(i));
%     end
% 
% 
%    fprintf(' %s \n','');
%    fprintf(' %s \n','');
%    fprintf(' %s \n','');
%    fprintf(' %s \n','-------------------');
%    fprintf(' %s \n','Bad Data Analysis');
%    fprintf(' %s \n','-----------------------------------------------------------------------');
%    fprintf(' %s \n','Measurement  Measurement    Caculated  Difference   Covariance   y_norm ');
%    fprintf(' %s \n','Number       Value          Value                   of error            ');
%    fprintf(' %s \n','-----------------------------------------------------------------------');
%    fprintf(' %s\n','Pflow measurements:');
%    for i=1:numFmeas
%        fprintf('%5d       %8.3f       %8.3f     %8.3f     %8.3f    %8.3f\n', i, FmeasPvalue(i)*baseMVA,  FcalcPvalue(i)*baseMVA, (FmeasPvalue(i)-FcalcPvalue(i))*baseMVA,sigmaPflow(i),y1_norm(i)*baseMVA);
%    end
%    fprintf(' %s\n','Qflow measurements:');
%    for i=1:numFmeas
%        fprintf('%5d       %8.3f       %8.3f     %8.3f     %8.3f    %8.3f\n', i, FmeasQvalue(i)*baseMVA,  FcalcQvalue(i)*baseMVA, (FmeasQvalue(i)-FcalcQvalue(i))*baseMVA,sigmaQflow(i),y2_norm(i)*baseMVA);
%    end 
%    fprintf(' %s\n','Pinjection measurements:');
%    for i=1:numImeas
%        fprintf('%5d       %8.3f       %8.3f     %8.3f     %8.3f    %8.3f\n', Imeasbus(i), ImeasPvalue(i)*baseMVA,  IcalcPvalue(i)*baseMVA, (ImeasPvalue(i)-IcalcPvalue(i))*baseMVA,sigmaPinjec(i),y3_norm(i)*baseMVA);
%    end
%    fprintf(' %s\n','Qinjection measurements:');
%    for i=1:numImeas
%        fprintf('%5d       %8.3f       %8.3f     %8.3f     %8.3f    %8.3f\n', Imeasbus(i), ImeasQvalue(i)*baseMVA,  IcalcQvalue(i)*baseMVA, (ImeasQvalue(i)-IcalcQvalue(i))*baseMVA,sigmaQinjec(i),y4_norm(i)*baseMVA);
%    end 
% 
% %Caculate all the flow and injection values from the voltages magnitudes and angles obtained from the state estimator. 
% for i=1:numbus
%     vr(i)=z.vmag(i)*(cos(z.vang(i))+sqrt(-1)*sin(z.vang(i)));
% end
% 
% 
% for i = 1:numline
%     a=frombus(i);
%     b=tobus(i);
%     Iline = (vr(a) - vr(b))/(R(i)+sqrt(-1)*X(i));
%     Iline = Iline + vr(a)*(0.0 + sqrt(-1)*Bcap(i)*0.5);
%     Sline = vr(a) * conj(Iline);
%     FcalcPvalue(i) = real(Sline);
% end
% for i=1:numline
%     a=tobus(i);
%     b=frombus(i);
%     Iline = (vr(a) - vr(b))/(R(i)+sqrt(-1)*X(i));
%     Iline = Iline + vr(a)*(0.0 + sqrt(-1)*Bcap(i)*0.5);
%     Sline = vr(a) * conj(Iline);
%     FcalcPvalue(i+numline) = real(Sline);
% end
% for i = 1:numline
%     a=frombus(i);
%     b=tobus(i);
%     Iline = (vr(a) - vr(b))/(R(i)+sqrt(-1)*X(i));
%     Iline = Iline + vr(a)*(0.0 + sqrt(-1)*Bcap(i)*0.5);
%     Sline = vr(a) * conj(Iline);
%     FcalcQvalue(i)= imag(Sline);
% end
% 
%  for i = 1:numline
%     a=tobus(i);
%     b=frombus(i);
%     Iline = (vr(a) - vr(b))/(R(i)+sqrt(-1)*X(i));
%     Iline = Iline + vr(a)*(0.0 + sqrt(-1)*Bcap(i)*0.5);
%     Sline = vr(a) * conj(Iline);
%     FcalcQvalue(i+numline)= imag(Sline);
% end
% IcalcPvalue=zeros(numbus,1);
% IcalcQvalue=zeros(numbus,1);
% for i=1:numbus
%     for b=1:numbus
%         Pinj =z.vmag(i)*z.vmag(b)*(cos(z.vang(i)-z.vang(b))*G(i,b)-sin(z.vang(i)-z.vang(b))*B(i,b));
%         IcalcPvalue(i)=IcalcPvalue(i)+Pinj;
%     end
% end
% 
% for i=1:numbus
%     for b=1:numbus
%         Qinj = z.vmag(i)*z.vmag(b)*(sin(z.vang(i)-z.vang(b))*G(i,b)+cos(z.vang(i)-z.vang(b))*B(i,b));
%         IcalcQvalue(i)=IcalcQvalue(i)+Qinj;
%     end
% end
% 
% 
%     % Display AC state estimator result
%      fprintf(' %s \n',' ');
%      fprintf(' %s \n','State Estimator Solution:');
%      fprintf(' %s \n',' ');
%      fprintf(' %s \n','Bus       Volt.    Volt.    Mw        Mvar      To      Mw        Mvar       ');
%      fprintf(' %s \n','Number    Mag.     Angle    Injec.    Injec.    Bus     Flow      Flow       ');
%      fprintf(' %s \n','------    ------   ------  --------  --------   ----   -------   -------     ');          
% 
%     for i=1:numbus
%         fprintf('%5d  %8.3f  %8.3f  %8.3f  %8.3f \n',i, z.vmag(i),z.vang(i),IcalcPvalue(i)*baseMVA,IcalcQvalue(i)*baseMVA); 
%         for j=1:numline
%             a=frombus(j);
%             b=tobus(j);
%            if b==i
%               fprintf('                                               %5d    %8.3f  %8.3f   \n', a, FcalcPvalue(j+numline)*baseMVA, FcalcQvalue(j+numline)*baseMVA);
%            elseif a==i
%               fprintf('                                               %5d    %8.3f  %8.3f   \n', b, FcalcPvalue(j)*baseMVA, FcalcQvalue(j)*baseMVA);
%            end
%         end
%     end