%Edit Injection Measurement


      fprintf('%s \n', 'which Injection measurement data you would like to edit?');
      fprintf('%s \n', 'P Injection measurement data:');
      for i=1:numImeas
          fprintf('%5d %8.3f\n',Imeasbus(i), ImeasPvalue(i)*baseMVA);
      end
       fprintf('%s \n', 'Q Injection measurement data:');
      for i=1:numImeas
          fprintf('%5d %8.3f\n',Imeasbus(i), ImeasQvalue(i)*baseMVA);
      end
      Injectionbus=input('Injection Measurement Bus Number:','s');
      [x status]=str2num(Injectionbus);
      if status==0
          fprintf('%s\n','No change has been made');
          edit_measurement;  
      else
          Injectionbus=str2num(Injectionbus);
          for i=1:numImeas
              if Imeasbus(i,1)==Injectionbus;
                  inj_num=i;   % Injection bus index
              end
          end
          for i=1:numImeas
              if Injectionbus==Imeasbus(i)
                  break 
              else
                  i=i+1;
              end
          end
          while i>numImeas
                 fprintf('%s\n','ERROR!!!');
                 fprintf('%s\n','Please enter the bus number again!');
                 Injectionbus=input('Injection Measurement Bus Number:');
                 for i=1:numImeas
                     if Imeasbus(i,1)==Injectionbus;
                        inj_num=i;   % Injection bus index
                     end
                 end
                 for i=1:numImeas
                     if Injectionbus==Imeasbus(i)
                        break 
                     else
                        i=i+1;
                     end
                 end
          end
          fprintf('%s \n', 'Please input new values to the injection measurement data (MW) of the bus you have selected');
          IinjPnew=0;
          IinjQnew=0;
          IinjPnew= input(' IinjPnew: ','s');
          [x status]=str2num(IinjPnew);
          if status==0
              fprintf('%s\n','No change has been made to P injection measurement data');  
          else
              IinjPnew=str2num(IinjPnew);
              ImeasPvalue(inj_num)=IinjPnew/baseMVA; 
          end
          IinjQnew= input(' IinjQnew: ','s');
          [x status]=str2num(IinjQnew);
          if status==0
              fprintf('%s\n','No change has been made to Q injection measurement data');  
          else
              IinjQnew=str2num(IinjQnew);
              ImeasQvalue(inj_num)=IinjQnew/baseMVA; 
          end
          fprintf('%s \n','Now all the injection measurement become:');
          fprintf('%s \n', 'P Injection measurement data:');
          for i=1:numImeas
              fprintf('%5d %8.3f\n',Imeasbus(i), ImeasPvalue(i)*baseMVA);
          end
           fprintf('%s \n', 'Q Injection measurement data:');
          for i=1:numImeas
              fprintf('%5d %8.3f\n',Imeasbus(i), ImeasQvalue(i)*baseMVA);
          end
          fprintf(' %s \n',' ');

          edit_measurement;
      end