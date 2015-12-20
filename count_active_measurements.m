%num_active_measurements is the total number of measurements being
% used it will be <= num_measurements depending on how many measurements
% have status set to 0
num_active_measurements = 0;
% num_active_measurements
% count active Voltage Measurements

    for i_vmeas = 1:numVmeas
        if Vmeasbus(i_vmeas) ~= refbus
            if Vmeasstatus(i_vmeas) == 1
                 num_active_measurements = num_active_measurements + 1;
            end
        elseif Vmeasbus(i_vmeas) == refbus
             num_active_measurements = num_active_measurements + 1;
        end
%         i_vmeas
%         num_active_measurements
    end
% num_active_measurements
    % Count active Angle Measurements
    for i_ameas = 1:numAmeas
      if Ameasbus(i_ameas) ~= refbus  
        if Ameasstatus(i_ameas) == 1
            num_active_measurements = num_active_measurements + 1;
        end
      elseif Ameasbus(i_ameas) == refbus 
            num_active_measurements = num_active_measurements + 1;
      end
    end
% num_active_measurements
% count active Injection Measurements
    for i_injmeas=1:numImeas
        if Imeasbus(i_injmeas) ~= refbus
             if Imeasstatus(i_injmeas) == 1
                 num_active_measurements = num_active_measurements + 2;
            end;
        end
    end
% num_active_measurements
% Count active flow measurements
    for i_flowmeas = 1:numFmeas
        if Fmeasstatus(i_flowmeas) == 1
            num_active_measurements = num_active_measurements + 2;
        end
    end
%     num_active_measurements
% num_active_measurements