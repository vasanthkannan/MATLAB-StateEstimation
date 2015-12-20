measurement_index = 1;

    % Caculate Residual 
    meas_error = zeros(num_active_measurements, 1);

    J = 0.0;

    % Build the error vector rows for Voltage Measurements
    Vmeas_estimated_value = zeros(numVmeas,1);
    for i_vmeas = 1:numVmeas
        i = Vmeasbus(i_vmeas);
        Vmeas_estimated_value(i_vmeas) = Vmag(i);
        
        if i ~= refbus
            
            if Vmeasstatus(i_vmeas) == 1
                meas_error(measurement_index)= Vmeas_measurement_value(i_vmeas) - Vmeas_estimated_value(i_vmeas);
                J = J + (Vmeas_measurement_value(i_vmeas)-Vmeas_estimated_value(i_vmeas))^2/Vmeassigma(i_vmeas)^2;
                measurement_index = measurement_index + 1;
            end
            
        elseif i == refbus
                   
                Vmeas_estimated_value(i_vmeas) = Vmeas_measurement_value(i_vmeas);
                meas_error(measurement_index) = 0.0;
                measurement_index = measurement_index + 1;
   
        end
 
    end
    


    % Build the error vector rows for Angle Measurements
    Ameas_estimated_value = zeros(numAmeas,1);
    for i_ameas = 1:numAmeas
        i = Ameasbus(i_ameas);
        Ameas_estimated_value(i_ameas) = Theta(i);
        
        if i ~= refbus
        
            if Ameasstatus(i_ameas) == 1
                meas_error(measurement_index)= Ameas_measurement_value(i_ameas) - Ameas_estimated_value(i_ameas);
                J = J + (Ameas_measurement_value(i_ameas)-Ameas_estimated_value(i_ameas))^2/Ameassigma(i_ameas)^2;
                measurement_index = measurement_index + 1;
            end
            
        elseif i == refbus
     
                Ameas_estimated_value(i_ameas) = Ameas_measurement_value(i_ameas);
                meas_error(measurement_index)= 0.0;
                measurement_index = measurement_index + 1;
        end     
        
    end


    % Residual for Injection Measurements
    Imeas_estimated_Pvalue=zeros(numImeas,1);
    Imeas_estimated_Qvalue=zeros(numImeas,1);
    for i_injmeas=1:numImeas
        
         i = Imeasbus(i_injmeas);
         
         if i ~=refbus
         
             Pinj(i)=0;
             Qinj(i)=0;
             for j = 1:numbus
                if i == j
                 Pinj(i) = Pinj(i) + G(i,i)*Vmag(i)^2;
                 Qinj(i) = Qinj(i) - B(i,i)*Vmag(i)^2;
                else         
                 Pinj(i) = Pinj(i) + Vmag(i)*Vmag(j)*(G(i,j)*cos(Theta(i)-Theta(j))+B(i,j)*sin(Theta(i)-Theta(j))); 
                 Qinj(i) = Qinj(i) + Vmag(i)*Vmag(j)*(G(i,j)*sin(Theta(i)-Theta(j))-B(i,j)*cos(Theta(i)-Theta(j)));
                end
             end
             Imeas_estimated_Pvalue(i_injmeas)= Pinj(i);
             Imeas_estimated_Qvalue(i_injmeas)= Qinj(i);
             if Imeasstatus(i_injmeas) == 1

                 J = J + (Imeas_measurement_Pvalue(i_injmeas)-Imeas_estimated_Pvalue(i_injmeas))^2/ImeasPsigma(i_injmeas)^2;
                 meas_error(measurement_index)= Imeas_measurement_Pvalue(i_injmeas)-Imeas_estimated_Pvalue(i_injmeas);
                 measurement_index = measurement_index + 1;


                 J = J + (Imeas_measurement_Qvalue(i_injmeas)-Imeas_estimated_Qvalue(i_injmeas))^2/ImeasQsigma(i_injmeas)^2;
                 meas_error(measurement_index)= Imeas_measurement_Qvalue(i_injmeas)-Imeas_estimated_Qvalue(i_injmeas);
                 measurement_index = measurement_index + 1;

            end;
         end
    end


    % Residual for flow measurements
    Fmeas_estimated_Pvalue = zeros(numFmeas,1);
    Fmeas_estimated_Qvalue = zeros(numFmeas,1);
    for i_flowmeas = 1:numFmeas
        i=Fmeasfrombus(i_flowmeas);
        j=Fmeastobus(i_flowmeas);
        ibranch = Fmeasbranch(i_flowmeas);
        Gline =  R(ibranch)/(R(ibranch)^2 + X(ibranch)^2);
        Bline = -X(ibranch)/(R(ibranch)^2 + X(ibranch)^2);
        Bcapline = Bcap(ibranch)/2;

        Fmeas_estimated_Pvalue(i_flowmeas) = Gline*Vmag(i)^2 - Gline*Vmag(i)*Vmag(j)*cos(Theta(i)-Theta(j))-Bline*Vmag(i)*Vmag(j)*sin(Theta(i)-Theta(j));
        Fmeas_estimated_Qvalue(i_flowmeas)= -1*Vmag(i)*Vmag(j)*sin(Theta(i)-Theta(j))*Gline - Vmag(i)^2*Bline + Vmag(i)*Vmag(j)*cos(Theta(i)-Theta(j))*Bline - (Vmag(i)^2)*Bcapline;

        if Fmeasstatus(i_flowmeas) == 1
            J = J + (Fmeas_measurement_Pvalue(i_flowmeas)- Fmeas_estimated_Pvalue(i_flowmeas))^2/FmeasPsigma(i_flowmeas)^2;
            meas_error(measurement_index)= Fmeas_measurement_Pvalue(i_flowmeas)- Fmeas_estimated_Pvalue(i_flowmeas);
            measurement_index = measurement_index + 1;

            J = J + (Fmeas_measurement_Qvalue(i_flowmeas)-Fmeas_estimated_Qvalue(i_flowmeas))^2/FmeasQsigma(i_flowmeas)^2;
            meas_error(measurement_index)= Fmeas_measurement_Qvalue(i_flowmeas)-Fmeas_estimated_Qvalue(i_flowmeas);
            measurement_index = measurement_index + 1;

        end
    end
   