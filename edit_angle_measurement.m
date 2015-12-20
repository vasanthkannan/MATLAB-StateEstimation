% edit_angle_measurement

      fprintf('%s \n', 'which angle measurement would you like to edit?');
      fprintf('%s    %s     %s \n', 'angle meas number','Ameasbus','Ameasvalue');
      for i=1:numAmeas
          fprintf('%5d             %5d        %8.3f \n', i,Ameasbus(i),Ameasvalue(i));
      end
       anglenumber=input('angle measurement number:','s');
      [x status]=str2num(anglenumber);
      if status==0
          fprintf('%s\n','No change has been made');
          edit_measurement;  
      else
          anglenumber=str2num(anglenumber);
          while (anglenumber>numAmeas||anglenumber<=0)
                 fprintf('%s\n','ERROR!!!');
                 fprintf('%s \n', 'which angle measurement would you like to edit?');
                 anglenumber=input('angle measurement number:');
          end
          a_bus=Ameasbus(anglenumber);
          fprintf('%s %5d \n', 'Volt bus:',a_bus);
          % Later I should add error signal if the bus number is not right.
          fprintf('%s \n', 'Please input a new value to the angle measurement data for the bus you have selected');
          Ameasnew = input(' Ameasnew: ','s');
          [x status]=str2num(Ameasnew);
          if status==0
              fprintf('%s\n','No change has been made to Voltage measurement data');  
          else
              Ameasnew=str2num(Ameasnew);
              Ameasvalue(anglenumber)=Ameasnew; 
          end
          
          fprintf('%s \n','Now the new angle data become:');
          fprintf('%s    %s     %s\n', 'angle meas number','Ameasbus','Ameasvalue');
          for i=1:numAmeas
              fprintf('%5d               %5d        %8.3f \n', i,Ameasbus(i),Ameasvalue(i));
          end
          edit_measurement;
      end