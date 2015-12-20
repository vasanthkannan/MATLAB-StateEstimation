% edit_flow_measurement

      fprintf('%s \n', 'which flow measurement data you would like to edit?');
      fprintf('%s   %s  %s    %s    %s     %s\n', 'flow number','branch number','Fmeasfrombus','Fmeastobus','FmeasPflow','FmeasQflow');
      for i=1:numFmeas
          fprintf('%5d        %5d            %5d            %5d        %8.3f        %8.3f \n', i,Fmeasbranch(i),Fmeasfrombus(i),Fmeastobus(i),FmeasPvalue(i)*baseMVA,FmeasQvalue(i)*baseMVA);
      end
       flownumber=input('flow measurement number:','s');
      [x status]=str2num(flownumber);
      if status==0
          fprintf('%s\n','No change has been made');
          edit_measurement;  
      else
          flownumber=str2num(flownumber);
          while (flownumber>numFmeas||flownumber<=0)
                 fprintf('%s\n','ERROR!!!');
                 fprintf('%s \n', 'which flow measurement data you would like to edit?');
                 flownumber=input('flow measurement number:');
          end
          a=Fmeasfrombus(flownumber);
          b=Fmeastobus(flownumber);
          fprintf('%s %5d    %s %5d \n', 'From bus:',a,'To bus:',b);
          % 
          fprintf('%s \n', 'Please input new values to the flow measurement data of the line you have selected');
          FmeasPnew= input(' FmeasPnew: ','s');
          [x status]=str2num(FmeasPnew);
          if status==0
              fprintf('%s\n','No change has been made to P flow measurement data');  
          else
              FmeasPnew=str2num(FmeasPnew);
              FmeasPvalue(flownumber)=FmeasPnew/baseMVA; 
          end
          FmeasQnew= input(' FmeasQnew: ','s'); 
          [x status]=str2num(FmeasQnew);
          if status==0
              fprintf('%s\n','No change has been made to Q flow measurement data');  
          else
              FmeasQnew=str2num(FmeasQnew);
              FmeasQvalue(flownumber)=FmeasQnew/baseMVA; 
          end
          fprintf('%s \n','Now the new flow data become:');
          fprintf('%s   %s   %s   %s    %s     %s\n', 'flow number','branch number','Fmeasfrombus','Fmeastobus','FmeasPflow','FmeasQflow');
          for i=1:numFmeas
              fprintf('%5d        %5d            %5d           %5d        %8.3f        %8.3f \n', i,Fmeasbranch(i),Fmeasfrombus(i),Fmeastobus(i),FmeasPvalue(i)*baseMVA,FmeasQvalue(i)*baseMVA);
          end
          edit_measurement;
      end