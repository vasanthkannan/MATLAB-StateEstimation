%scada

% Add measurement value tables to measurement data structures to 
% contain the measurement values

SCADA_voltage_meas_data.voltage_meas_value =  Vmeasvalue;
SCADA_voltage_meas_data.voltage_meas_status = Vmeasstatus;

SCADA_angle_meas_data.angle_meas_value =  Ameasvalue;
SCADA_angle_meas_data.angle_meas_status = Ameasstatus;

SCADA_injection_meas_data.injection_meas_Pvalue = ImeasPvalue;
SCADA_injection_meas_data.injection_meas_Qvalue = ImeasQvalue;
SCADA_injection_meas_data.injection_meas_status = Imeasstatus;

SCADA_flow_meas_data.flow_meas_Pvalue = FmeasPvalue;
SCADA_flow_meas_data.flow_meas_Qvalue = FmeasQvalue;
SCADA_flow_meas_data.flow_meas_status = Fmeasstatus;

%  Save StateEstimatorData.mat file to send to State Estimator

save('StateEstimatorData','SCADA_voltage_meas_data',...
                          'SCADA_angle_meas_data',...
                          'SCADA_injection_meas_data',...
                          'SCADA_flow_meas_data');who -file StateEstimatorData

