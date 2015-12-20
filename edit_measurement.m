%edit_measurement

k = menu('Edit Measurement Data','1) Edit Voltage Measurement',...
                                 '2) Edit Angle Measurement',...
                                 '3) Edit Injection Measurement',...
                                 '4) Edit Flow Measurement',...
                                 '5) Exit Editor'); 



if k==1
    edit_volt_measurement;
elseif k==2
    edit_angle_measurement;
elseif k==3
    edit_injection_measurement;
elseif k==4
    edit_flow_measurement;
elseif k==5
    Scadamenu;
end