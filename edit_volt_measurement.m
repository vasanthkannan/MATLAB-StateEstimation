% edit_volt_measurement

      fprintf('%s \n', 'which volt measurement would you like to edit?');
      fprintf('%s    %s     %s \n',...
              'volt meas number','Vmeasbus','Vmeasvalue');
      for i=1:numVmeas
          fprintf('%5d             %5d          %8.3f \n', i,Vmeasbus(i),Vmeasvalue(i));
      end
       voltnumber=input('volt measurement number:','s');
      [x status]=str2num(voltnumber);
      if status==0
          fprintf('%s\n','No change has been made');
          edit_measurement;  
      else
          voltnumber=str2num(voltnumber);
          while (voltnumber>numVmeas||voltnumber<=0)
                 fprintf('%s\n','ERROR!!!');
                 fprintf('%s \n', 'which volt measurement would you like to edit?');
                 voltnumber=input('volt measurement number:');
          end
          v_bus=Vmeasbus(voltnumber);
          fprintf('%s %5d \n', 'Volt bus:',v_bus);
          % Later I should add error signal if the bus number is not right.
          fprintf('%s \n', 'Please input a new value to the voltage measurement data for the bus you have selected');
          Vmeasnew = input(' Vmeasnew: ','s');
          [x status]=str2num(Vmeasnew);
          if status==0
              fprintf('%s\n','No change has been made to Voltage measurement data');  
          else
              Vmeasnew=str2num(Vmeasnew);
              Vmeasvalue(voltnumber)=Vmeasnew; 
          end
          
          fprintf('%s \n','Now the new voltage data become:');
          fprintf('%s    %s     %s\n', 'volt meas number','Vmeasbus','Vmeasvalue');
          for i=1:numVmeas
              fprintf('%5d             %5d          %8.3f \n', i,Vmeasbus(i),Vmeasvalue(i));
          end
          edit_measurement;
      end