measurement_index = 1;

% Voltage Measurements
for i_vmeas = 1:numVmeas
    i = Vmeasbus(i_vmeas);

    if Vmeasstatus(i_vmeas) == 1

        if i == Slack

        end    
    end
    measurement_index = measurement_index + 1;
end


% Angle Measurements
for i_ameas = 1:numAmeas
    i = Ameasbus(i_ameas);

    if Ameasstatus(i_ameas) == 1

        if i == Slack

        end
    end
    measurement_index = measurement_index + 1;
end


% Injection Measurements

for i_injmeas=1:numImeas
     i = Imeasbus(i_injmeas);

     if Imeasstatus(i_injmeas) == 1
         measurement_index = measurement_index + 1;



         measurement_index = measurement_index + 1;
    else
         measurement_index = measurement_index + 2
    end;
end


% Flow Measurements

for i_flowmeas = 1:numFmeas
    i=Fmeasfrombus(i_flowmeas);
    j=Fmeastobus(i_flowmeas);
    ibranch = Fmeasbranch(i_flowmeas);


    if Fmeasstatus(i_flowmeas) == 1

        measurement_index = measurement_index + 1;


        measurement_index = measurement_index + 1;
    else
        measurement_index = measurement_index + 2;
    end
end